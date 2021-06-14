use crate::{BoundingBox, Object, PrimitiveParameters, RealField};
use num_traits::Float;

/// Twister will twist an object by rotating it along the Z-Axis.
#[derive(Clone, Debug)]
pub struct Twister<S: RealField> {
    object: Box<dyn Object<S>>,
    height_scaler: S, // 2 * pi / (height for full rotation)
    value_scaler: S,
    bbox: BoundingBox<S>,
}

impl<S: RealField + From<f32> + Float + ::num_traits::FloatConst> Object<S> for Twister<S> {
    fn approx_value(&self, p: &na::Point3<S>, slack: S) -> S {
        let approx = self.bbox.distance(p);
        if approx <= slack {
            self.object
                .approx_value(&self.twist_point(p), slack / self.value_scaler)
                * self.value_scaler
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
        self.untwist_normal(&self.object.normal(&self.twist_point(&p)), p)
    }
}

impl<S: RealField + Float + ::num_traits::FloatConst + From<f32>> Twister<S> {
    /// Create a twisted version ob o.
    /// o: Object to be twisted, h: height for one full rotation
    pub fn new(o: Box<dyn Object<S>>, h: S) -> Self {
        let _2pi: S = S::PI() * From::from(2.);
        let mx = Float::max(Float::abs(o.bbox().min.x), Float::abs(o.bbox().max.x));
        let my = Float::max(Float::abs(o.bbox().min.y), Float::abs(o.bbox().max.y));
        let r = Float::hypot(mx, my);

        // The ratio of height and circumference (slope on the outer edge).
        let tan_a = Float::abs(h) / (_2pi * r);
        // The scaler is 1 / sin(a)
        // sin(atan(x)) =   x / sqrt(x^2 + 1)
        let one: S = From::from(1f32);
        let scaler = tan_a / Float::sqrt(tan_a * tan_a + one);

        let bbox = BoundingBox::<S>::new(
            &na::Point3::new(-r, -r, o.bbox().min.z),
            &na::Point3::new(r, r, o.bbox().max.z),
        );
        Twister {
            object: o,
            height_scaler: _2pi / h,
            value_scaler: scaler,
            bbox,
        }
    }
    fn twist_point(&self, p: &na::Point3<S>) -> na::Point3<S> {
        let p2 = ::na::Point2::new(p.x, p.y);
        let angle = p.z * self.height_scaler;
        type Rota<S> = ::na::Rotation<S, 2>;
        let trans = Rota::new(angle);
        let rp2 = trans.transform_point(&p2);
        na::Point3::new(rp2.x, rp2.y, p.z)
    }
    // Apply tilt to the vector.
    // Since Surfaces are twisted, all normals will be tilted, depending on the radius.
    fn tilt_normal(&self, normal: na::Vector3<S>, p: &na::Point3<S>) -> na::Vector3<S> {
        let radius_v = ::na::Vector2::new(p.x, p.y);
        let radius = radius_v.norm();
        let radius_v = radius_v / radius;
        // Calculate tangential unit na::Vector3<S> at p.
        let tangent_v = ::na::Vector2::new(radius_v.y, -radius_v.x);

        // Project the in plane component of normal onto tangent.
        let planar_normal = ::na::Vector2::new(normal.x, normal.y);
        let tangential_projection = tangent_v.dot(&planar_normal);

        // Calculate the shear at p.
        let tangential_shear = radius * self.height_scaler;

        // Subtract from normal.z.
        let mut result = normal;
        result.z -= tangential_shear * tangential_projection;

        // Normalize.
        result.normalize()
    }
    fn untwist_normal(&self, v: &na::Vector3<S>, p: &na::Point3<S>) -> na::Vector3<S> {
        let v2 = ::na::Vector2::new(v.x, v.y);
        let angle = -p.z * self.height_scaler;
        let trans = ::na::Rotation2::new(angle);
        let rv2 = trans.transform_vector(&v2);
        self.tilt_normal(na::Vector3::new(rv2.x, rv2.y, v.z), p)
    }
}

#[cfg(test)]
mod test {
    use crate::test::MockObject;
    use crate::*;

    #[test]
    fn simple() {
        let m = MockObject::new_with_bbox(
            10.0,
            na::Vector3::new(1., 0., 0.),
            BoundingBox::new(
                &na::Point3::new(-1., -1., -100.),
                &na::Point3::new(1., 1., 100.),
            ),
        );
        let t = Twister::new(Box::new(m), 4.);
        assert_relative_eq!(
            t.approx_value(&na::Point3::new(0., 0., 0.), 0.),
            4.104_846_065_998_354
        );
        assert_relative_eq!(
            t.normal(&na::Point3::new(1., 0., 0.)),
            na::Vector3::new(1., 0., 0.)
        );

        assert_relative_eq!(
            t.approx_value(&na::Point3::new(0., 0., 1.), 0.),
            4.104_846_065_998_354
        );
        ulps_eq!(
            t.normal(&na::Point3::new(1., 0., 1.)),
            na::Vector3::new(0., -0.537_029_272_146_315_1, 0.843_563_608_068_768_6)
        );

        assert_relative_eq!(
            t.approx_value(&na::Point3::new(0., 0., 2.), 0.),
            4.104_846_065_998_354
        );
        assert_relative_eq!(
            t.normal(&na::Point3::new(1., 0., 2.)),
            na::Vector3::new(-1., 0., 0.)
        );

        assert_relative_eq!(
            t.approx_value(&na::Point3::new(0., 0., 3.), 0.),
            4.104_846_065_998_354
        );
        assert_relative_eq!(
            t.normal(&na::Point3::new(1., 0., 3.)),
            na::Vector3::new(0., 0.537_029_272_146_315_1, 0.843_563_608_068_768_6)
        );

        assert_relative_eq!(
            t.approx_value(&na::Point3::new(0., 0., 4.), 0.),
            4.104_846_065_998_354
        );
        assert_relative_eq!(
            t.normal(&na::Point3::new(1., 0., 4.)),
            na::Vector3::new(1., 0.0, 0.0),
            epsilon = 10e-10
        );
    }
}
