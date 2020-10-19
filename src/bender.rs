use crate::{BoundingBox, Object, PrimitiveParameters, RealField};
use num_traits::{Float, FloatConst};

/// Bender create an implicit function that represents a bended version of it's input.
/// The object will be bend around the Z-Axis.
/// E.g. bending a cylinder along the X-Axis (and translated away from the Z-Axis) will result in a Torus.
#[derive(Clone, Debug)]
pub struct Bender<S: RealField> {
    object: Box<dyn Object<S>>,
    width_scaler: S, // width_for_full_rotation / (2. * PI),
    bbox: BoundingBox<S>,
}

impl<S: RealField + From<f32> + Float + ::num_traits::FloatConst> Object<S> for Bender<S> {
    fn approx_value(&self, p: &na::Point3<S>, slack: S) -> S {
        let approx = self.bbox.distance(p);
        if approx <= slack {
            let mut obj_p = self.to_polar(p);
            let r = obj_p.y;

            // If the bended object is a ring, and p is in the center, return the distance to inner
            // margin (bbox.min.y) of the (bent) bounding box.
            let center_to_bbox = self.object.bbox().min.y - r;
            if center_to_bbox > slack {
                return center_to_bbox;
            }

            // let circumference = 2. * PI * r;
            // let width_for_full_rotation = self.width_scaler * 2. * PI;
            // let x_scale = circumference / width_for_full_rotation;
            let x_scale = r / self.width_scaler;
            let x_scaler = Float::min(x_scale, From::from(1f32));

            obj_p.x *= self.width_scaler;
            self.object.approx_value(&obj_p, slack / x_scaler) * x_scaler
        } else {
            approx
        }
    }
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
    fn set_parameters(&mut self, p: &PrimitiveParameters<S>) {
        self.object.set_parameters(p);
    }
    fn normal(&self, p: &na::Point3<S>) -> na::Vector3<S> {
        let polar_p = self.to_polar(p);
        let mut obj_p = polar_p;
        obj_p.x *= self.width_scaler;
        self.bend_normal(self.object.normal(&obj_p), polar_p)
    }
}

impl<S: RealField + Float + FloatConst + From<f32>> Bender<S> {
    /// Create a new bent object.
    /// o: Object to be bent, w: width (x) for one full rotation
    pub fn new(o: Box<dyn Object<S>>, w: S) -> Self {
        let bbox = BoundingBox::new(
            &na::Point3::new(-o.bbox().max.y, -o.bbox().max.y, o.bbox().min.z),
            &na::Point3::new(o.bbox().max.y, o.bbox().max.y, o.bbox().max.z),
        );
        let _2pi: S = S::PI() * From::from(2.);
        Bender {
            object: o,
            width_scaler: w / _2pi,
            bbox,
        }
    }
    fn to_polar(&self, p: &na::Point3<S>) -> na::Point3<S> {
        let phi = Float::atan2(p.x, -p.y);
        let r = Float::hypot(p.x, p.y);
        na::Point3::new(phi, r, p.z)
    }
    fn tilt_normal(&self, mut normal: na::Vector3<S>, polar_p: na::Point3<S>) -> na::Vector3<S> {
        let r = polar_p.y;
        let _2pi: S = S::PI() * From::from(2.);
        let circumference = _2pi * r;
        let width_for_one_full_rotation = self.width_scaler * _2pi;
        let scale_along_x = circumference / width_for_one_full_rotation;
        normal.x /= scale_along_x;
        normal.normalize()
    }
    fn bend_normal(&self, v: na::Vector3<S>, polar_p: na::Point3<S>) -> na::Vector3<S> {
        let v = self.tilt_normal(v, polar_p);
        let phi = polar_p.x;
        let v2 = ::na::Vector2::new(v.x, -v.y);
        let trans = ::na::Rotation2::new(phi);
        let rv2 = trans.transform_vector(&v2);
        na::Vector3::new(rv2.x, rv2.y, v.z)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test::MockObject;

    #[test]
    fn values_in_quadrants() {
        let m = MockObject::new(10.0, na::Vector3::new(1., 0., 0.));
        let b = Bender::new(Box::new(m), 4.);

        assert_relative_eq!(b.approx_value(&na::Point3::new(0., 1., 0.), 0.), 10.);
        assert_relative_eq!(b.approx_value(&na::Point3::new(-1., 0., 0.), 0.), 10.);
        assert_relative_eq!(b.approx_value(&na::Point3::new(0., -1., 0.), 0.), 10.);
        assert_relative_eq!(b.approx_value(&na::Point3::new(1., 0., 0.), 0.), 10.);
    }
    #[test]
    fn normal_x_in_quadrants() {
        let m = MockObject::new(10.0, na::Vector3::new(1., 0., 0.));
        let b = Bender::new(Box::new(m), 4.);

        assert_relative_eq!(
            b.normal(&na::Point3::new(0., 1., 0.)),
            na::Vector3::new(-1., 0.0, 0.0)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(-1., 0., 0.)),
            na::Vector3::new(0., -1., 0.)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(0., -1., 0.)),
            na::Vector3::new(1., 0., 0.)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(1., 0., 0.)),
            na::Vector3::new(0., 1., 0.)
        );
    }
    #[test]
    fn normal_y_in_quadrants() {
        let m = MockObject::new(10.0, na::Vector3::new(0., 1., 0.));
        let b = Bender::new(Box::new(m), 4.);

        assert_relative_eq!(
            b.normal(&na::Point3::new(0., 1., 0.)),
            na::Vector3::new(0., 1., 0.0)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(-1., 0., 0.)),
            na::Vector3::new(-1., 0., 0.)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(0., -1., 0.)),
            na::Vector3::new(0., -1., 0.)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(1., 0., 0.)),
            na::Vector3::new(1., -0.0, 0.)
        );
    }
    #[test]
    fn normal_z_in_quadrants() {
        let m = MockObject::new(10.0, na::Vector3::new(0., 0., 1.));
        let b = Bender::new(Box::new(m), 4.);

        assert_relative_eq!(
            b.normal(&na::Point3::new(0., 1., 0.)),
            na::Vector3::new(0., 0., 1.)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(-1., 0., 0.)),
            na::Vector3::new(0., 0., 1.)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(0., -1., 0.)),
            na::Vector3::new(0., 0., 1.)
        );

        assert_relative_eq!(
            b.normal(&na::Point3::new(1., 0., 0.)),
            na::Vector3::new(0., 0., 1.)
        );
    }
}
