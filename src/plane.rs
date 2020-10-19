use crate::{BoundingBox, Object, RealField};
use num_traits::Float;

pub trait Axis: ::std::fmt::Debug + Clone + ::std::marker::Sync + ::std::marker::Send {
    fn value() -> usize;
    fn inverted() -> bool;
}

#[derive(Clone, Debug)]
pub struct AxisX {}
impl Axis for AxisX {
    fn value() -> usize {
        0
    }
    fn inverted() -> bool {
        false
    }
}

#[derive(Clone, Debug)]
pub struct AxisNegX {}
impl Axis for AxisNegX {
    fn value() -> usize {
        0
    }
    fn inverted() -> bool {
        true
    }
}

#[derive(Clone, Debug)]
pub struct AxisY {}
impl Axis for AxisY {
    fn value() -> usize {
        1
    }
    fn inverted() -> bool {
        false
    }
}

#[derive(Clone, Debug)]
pub struct AxisNegY {}
impl Axis for AxisNegY {
    fn value() -> usize {
        1
    }
    fn inverted() -> bool {
        true
    }
}

#[derive(Clone, Debug)]
pub struct AxisZ {}
impl Axis for AxisZ {
    fn value() -> usize {
        2
    }
    fn inverted() -> bool {
        false
    }
}

#[derive(Clone, Debug)]
pub struct AxisNegZ {}
impl Axis for AxisNegZ {
    fn value() -> usize {
        2
    }
    fn inverted() -> bool {
        true
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Plane<A: Axis, S: RealField> {
    distance_from_zero: S,
    bbox: BoundingBox<S>,
    normal: na::Vector3<S>,
    _phantom: ::std::marker::PhantomData<A>,
}

impl<A: Axis, S: From<f32> + RealField + Float> Plane<A, S> {
    pub fn new(distance_from_zero: S) -> Self {
        let d = distance_from_zero;
        let mut p_neg = na::Point3::new(S::neg_infinity(), S::neg_infinity(), S::neg_infinity());
        let mut p_pos = na::Point3::new(S::infinity(), S::infinity(), S::infinity());
        let zero: S = From::from(0f32);
        let mut normal = na::Vector3::new(zero, zero, zero);

        if A::inverted() {
            p_neg[A::value()] = -d;
            normal[A::value()] = From::from(-1.);
        } else {
            p_pos[A::value()] = d;
            normal[A::value()] = From::from(1.);
        }

        Plane {
            distance_from_zero: d,
            bbox: BoundingBox::new(&p_neg, &p_pos),
            normal,
            _phantom: ::std::marker::PhantomData,
        }
    }
}

impl<A: 'static + Axis, S: Float + From<f32> + RealField> Object<S> for Plane<A, S> {
    fn approx_value(&self, p: &na::Point3<S>, _: S) -> S {
        let p: S = p[A::value()];
        let ap: S = Float::abs(p);
        if Float::is_sign_positive(p) != A::inverted() {
            ap - self.distance_from_zero
        } else {
            -ap - self.distance_from_zero
        }
    }
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
    fn normal(&self, _: &na::Point3<S>) -> na::Vector3<S> {
        self.normal
    }
}

/// A YZ-Plane.
pub type PlaneX<S> = Plane<AxisX, S>;
/// A negative YZ-Plane.
pub type PlaneNegX<S> = Plane<AxisNegX, S>;
/// A XZ-Plane
pub type PlaneY<S> = Plane<AxisY, S>;
/// A negative XZ-Plane
pub type PlaneNegY<S> = Plane<AxisNegY, S>;
/// A XY-Plane
pub type PlaneZ<S> = Plane<AxisZ, S>;
/// A negative XZ-Plane
pub type PlaneNegZ<S> = Plane<AxisNegZ, S>;

/// An arbitrary (not axis aligned) plane.
#[derive(Clone, Debug, PartialEq)]
pub struct NormalPlane<S: RealField> {
    bbox: BoundingBox<S>,
    normal: na::Vector3<S>,
    p: S,
}

impl<S: From<f32> + RealField + Float> NormalPlane<S> {
    /// Create a plane in hessian form.
    pub fn from_normal_and_p(normal: na::Vector3<S>, p: S) -> Self {
        NormalPlane {
            bbox: BoundingBox::infinity(),
            normal,
            p,
        }
    }
    /// Create a plane from 3 points.
    pub fn from_3_points(a: &na::Point3<S>, b: &na::Point3<S>, c: &na::Point3<S>) -> Self {
        let v1 = a - c;
        let v2 = b - c;
        let normal = v1.cross(&v2).normalize();
        let p = normal.dot(&a.coords);
        NormalPlane {
            bbox: BoundingBox::infinity(),
            normal,
            p,
        }
    }
}

impl<S: Float + From<f32> + RealField> Object<S> for NormalPlane<S> {
    fn approx_value(&self, x0: &na::Point3<S>, _: S) -> S {
        self.normal.dot(&x0.coords) - self.p
    }
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
    fn normal(&self, _: &na::Point3<S>) -> na::Vector3<S> {
        self.normal
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn simple() {
        let px = PlaneX::new(10.);
        assert_ulps_eq!(px.approx_value(&na::Point3::new(0., 0., 0.), 0.), -10.);
        assert_ulps_eq!(px.approx_value(&na::Point3::new(10., 0., 0.), 0.), 0.);
        assert_ulps_eq!(px.approx_value(&na::Point3::new(20., 0., 0.), 0.), 10.);
        assert_ulps_eq!(
            px.approx_value(&na::Point3::new(20., 1000., 1000.), 0.),
            10.
        );
    }

    #[test]
    fn hessian_x() {
        let px = NormalPlane::from_normal_and_p(na::Vector3::new(1., 0., 0.), 0.);
        assert_ulps_eq!(px.approx_value(&na::Point3::new(0., 0., 0.), 0.), 0.);
        assert_ulps_eq!(px.approx_value(&na::Point3::new(0., 10., 100.), 0.), 0.);
        assert_ulps_eq!(px.approx_value(&na::Point3::new(10., 0., 0.), 0.), 10.);
    }

    #[test]
    fn hessian_xyz() {
        let p = NormalPlane::from_normal_and_p(na::Vector3::new(1., 1., 1.).normalize(), 1.);
        assert_ulps_eq!(p.approx_value(&na::Point3::new(0., 0., 0.), 0.), -1.);
        assert_ulps_eq!(
            p.approx_value(&na::Point3::new(1., 1., 1.), 0.),
            Float::sqrt(3.) - 1.0
        );
        assert_ulps_eq!(
            p.approx_value(&na::Point3::new(2., 2., 2.), 0.),
            Float::sqrt(12.) - 1.0
        );
    }

    #[test]
    fn hessian_3points_x() {
        let p = NormalPlane::from_3_points(
            &na::Point3::new(10., 0., 0.),
            &na::Point3::new(10., 1., 0.),
            &na::Point3::new(10., 0., 1.),
        );
        assert_ulps_eq!(p.approx_value(&na::Point3::new(0., 0., 0.), 0.), -10.);
        assert_ulps_eq!(p.approx_value(&na::Point3::new(1., 1., 1.), 0.), -9.0);
        assert_ulps_eq!(p.approx_value(&na::Point3::new(100., 100., 100.), 0.), 90.);
    }

    #[test]
    fn hessian_3points() {
        let p = NormalPlane::from_3_points(
            &na::Point3::new(0., 1., 1.),
            &na::Point3::new(1., 0., 1.),
            &na::Point3::new(1., 1., 0.),
        );
        assert_ulps_eq!(
            p.approx_value(&na::Point3::new(0., 0., 0.), 0.),
            -Float::sqrt(4. / 3.)
        );
        assert_ulps_eq!(
            p.approx_value(&na::Point3::new(1., 1., 1.), 0.),
            Float::sqrt(3.) - Float::sqrt(4. / 3.)
        );
        assert_ulps_eq!(
            p.approx_value(&na::Point3::new(2., 2., 2.), 0.),
            Float::sqrt(12.) - Float::sqrt(4. / 3.)
        );
    }
}
