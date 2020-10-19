use crate::{BoundingBox, Object, RealField};
use num_traits::Float;

/// A cylinder along the Z-Axis
#[derive(Clone, Debug, PartialEq)]
pub struct Cylinder<S: RealField> {
    radius: S,
    bbox: BoundingBox<S>,
}

impl<S: RealField + Float> Cylinder<S> {
    /// Create a new infinite Cylinder (along the Z-Axis) of radius r.
    pub fn new(r: S) -> Self {
        Cylinder {
            radius: r,
            bbox: BoundingBox::new(
                &na::Point3::new(-r, -r, S::neg_infinity()),
                &na::Point3::new(r, r, S::infinity()),
            ),
        }
    }
}

impl<S: ::std::fmt::Debug + RealField + From<f32> + Float> Object<S> for Cylinder<S> {
    fn approx_value(&self, p: &na::Point3<S>, slack: S) -> S {
        let approx = self.bbox.distance(p);
        if approx <= slack {
            let zero: S = From::from(0f32);
            let pv = na::Vector3::new(p.x, p.y, zero);
            pv.norm() - self.radius
        } else {
            approx
        }
    }
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
    fn normal(&self, p: &na::Point3<S>) -> na::Vector3<S> {
        let zero: S = From::from(0f32);
        let pv = na::Vector3::new(p.x, p.y, zero);
        pv.normalize()
    }
}

/// A cone along the Z-Axis
#[derive(Clone, Debug, PartialEq)]
pub struct Cone<S: RealField> {
    slope: S,
    distance_multiplier: S,
    offset: S,            // Offset the singularity from Z-zero
    normal_multiplier: S, // muliplier for the normal caclulation
    bbox: BoundingBox<S>,
}

impl<S: RealField + Float + From<f32>> Cone<S> {
    /// Create a new infinite Cone (along the Z-Axis) for a given slope and and offset from origin.
    pub fn new(slope: S, offset: S) -> Self {
        let one: S = From::from(1f32);
        Cone {
            slope,
            distance_multiplier: one / Float::sqrt(slope * slope + one), // cos(atan(slope))
            offset,
            normal_multiplier: slope / Float::sqrt(slope * slope + one), // sin(atan(slope))
            bbox: BoundingBox::infinity(),
        }
    }
}

impl<S: ::std::fmt::Debug + RealField + From<f32> + Float> Object<S> for Cone<S> {
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
    fn set_bbox(&mut self, bbox: &BoundingBox<S>) {
        self.bbox = bbox.clone();
    }
    fn approx_value(&self, p: &na::Point3<S>, _: S) -> S {
        let radius = Float::abs(self.slope * (p.z + self.offset));
        let zero: S = From::from(0f32);
        let pv = na::Vector3::new(p.x, p.y, zero);
        (pv.norm() - radius) * self.distance_multiplier
    }
    fn normal(&self, p: &na::Point3<S>) -> na::Vector3<S> {
        let s = Float::signum(p.z + self.offset);
        let zero: S = From::from(0f32);
        let mut pv = na::Vector3::new(p.x, p.y, zero);
        pv.normalize_mut();
        pv *= self.distance_multiplier;
        pv.z = -s * self.normal_multiplier;
        pv
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn cylinder() {
        let cyl = Cylinder::new(1.0);
        assert_ulps_eq!(cyl.approx_value(&na::Point3::new(0., 0., 0.), 0.), -1.);
        assert_ulps_eq!(cyl.approx_value(&na::Point3::new(1., 0., 0.), 0.), 0.);
        assert_ulps_eq!(cyl.approx_value(&na::Point3::new(0., 1., 0.), 0.), 0.);
        assert_ulps_eq!(cyl.approx_value(&na::Point3::new(0., 10., 0.), 0.), 9.);
        assert_ulps_eq!(cyl.approx_value(&na::Point3::new(0., 10., 1000.), 0.), 9.);
    }

    #[test]
    fn cone() {
        let c = Cone::new(2., 10.);
        assert_ulps_eq!(c.approx_value(&na::Point3::new(0., 0., -10.), 0.), 0.);
        assert_ulps_eq!(c.approx_value(&na::Point3::new(2., 0., -11.), 0.), 0.);
    }
}
