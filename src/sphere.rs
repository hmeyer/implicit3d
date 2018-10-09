use alga::general::Real;
use na;
use num_traits::Float;
use {BoundingBox, Object};

/// A sphere. Simple.
#[derive(Clone, Debug, PartialEq)]
pub struct Sphere<S: Real> {
    radius: S,
    bbox: BoundingBox<S>,
}

impl<S: Real + Float> Sphere<S> {
    /// Create a new sphere of radius r.
    pub fn new(r: S) -> Box<Sphere<S>> {
        Box::new(Sphere {
            radius: r,
            bbox: BoundingBox::new(&na::Point3::new(-r, -r, -r), &na::Point3::new(r, r, r)),
        })
    }
}

impl<S: ::std::fmt::Debug + Real + Float + From<f32>> Object<S> for Sphere<S> {
    fn approx_value(&self, p: &na::Point3<S>, slack: S) -> S {
        let approx = self.bbox.distance(p);
        if approx <= slack {
            na::Vector3::new(p.x, p.y, p.z).norm() - self.radius
        } else {
            approx
        }
    }
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
    fn normal(&self, p: &na::Point3<S>) -> na::Vector3<S> {
        na::Vector3::new(p.x, p.y, p.z).normalize()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn simple() {
        let s = Sphere::new(1.);
        assert_ulps_eq!(s.approx_value(&na::Point3::new(0., 0., 0.), 0.), -1.);
        assert_ulps_eq!(s.approx_value(&na::Point3::new(1., 0., 0.), 0.), 0.);
        assert_ulps_eq!(s.approx_value(&na::Point3::new(2., 0., 0.), 0.), 1.);
    }
}
