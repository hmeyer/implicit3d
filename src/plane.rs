use {BoundingBox, Object};
use alga::general::Real;
use na;
use num_traits::Float;

pub trait Axis
    : ::std::fmt::Debug + Clone + ::std::marker::Sync + ::std::marker::Send {
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
pub struct Plane<A: Axis, S: Real> {
    distance_from_zero: S,
    bbox: BoundingBox<S>,
    normal: na::Vector3<S>,
    _phantom: ::std::marker::PhantomData<A>,
}

impl<A: Axis, S: From<f32> + Real + Float> Plane<A, S> {
    pub fn new(distance_from_zero: S) -> Box<Plane<A, S>> {
        let d = distance_from_zero;
        let mut p_neg = na::Point3::new(S::neg_infinity(), S::neg_infinity(), S::neg_infinity());
        let mut p_pos = na::Point3::new(S::infinity(), S::infinity(), S::infinity());
        let _0: S = From::from(0f32);
        let mut normal = na::Vector3::new(_0, _0, _0);


        if A::inverted() {
            p_neg[A::value()] = -d;
            normal[A::value()] = From::from(-1.);
        } else {
            p_pos[A::value()] = d;
            normal[A::value()] = From::from(1.);
        }


        Box::new(Plane {
            distance_from_zero: d,
            bbox: BoundingBox::new(&p_neg, &p_pos),
            normal: normal,
            _phantom: ::std::marker::PhantomData,
        })
    }
}

impl<A: 'static + Axis, S: Float + From<f32> + Real> Object<S> for Plane<A, S> {
    fn approx_value(&self, p: &na::Point3<S>, _: S) -> S {
        let p: S = p[A::value()];
        let ap: S = Float::abs(p);
        if Float::is_sign_positive(p) != A::inverted() {
            return ap - self.distance_from_zero;
        } else {
            return -ap - self.distance_from_zero;
        }
    }
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
    fn normal(&self, _: &na::Point3<S>) -> na::Vector3<S> {
        return self.normal.clone();
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
