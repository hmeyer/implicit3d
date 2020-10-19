//! ```implicit3d``` is a crate for creating
//! [3d implicit functions](https://en.wikipedia.org/wiki/Implicit_function).
//! Implicit functions evaluate to a scalar value for each point the 3d space.
//! They can be used to described object surfaces. If the function evaluates to negative values
//! the point is in the object, if the function evaluates positve this is outside the object.
//! If the function evaluates to zero the point is on the object surface.
//! This library allows to create implicit functions for 3d primitives (sphere, cylinder, cone,
//! box). Those primitives can be combined using
//! [CSG](https://en.wikipedia.org/wiki/Constructive_solid_geometry) and transformed.
//!
//! # Examples
//!
//! Create a Sphere:
//!
//! ```rust,no_run
//! let sphere = implicit3d::Sphere::new(1.0);
//! ```
//! Create a rounded Cube (as rounded intersection of 6 planes):
//!
//! ```rust,no_run
//! use std::fs::OpenOptions;
//! let px = Box::new(implicit3d::PlaneX::new(1.0));
//! let pnx = Box::new(implicit3d::PlaneNegX::new(1.0));
//! let py = Box::new(implicit3d::PlaneY::new(1.0));
//! let pny = Box::new(implicit3d::PlaneNegY::new(1.0));
//! let pz = Box::new(implicit3d::PlaneZ::new(1.0));
//! let pnz = Box::new(implicit3d::PlaneNegZ::new(1.0));
//! let cube = implicit3d::Intersection::from_vec(vec![px, pnx, py, pny, pz, pnz], 0.2);
//! ```

#![warn(missing_docs)]

#[cfg(test)]
#[macro_use]
extern crate approx;

extern crate nalgebra as na;

pub use bbox::BoundingBox;
use num_traits::Float;
use std::fmt::Debug;

/// A Combination of alga::general::RealField and na::RealField.
pub trait RealField: alga::general::RealField + na::RealField {}
impl RealField for f64 {}
impl RealField for f32 {}

mod transformer;
pub use self::transformer::AffineTransformer;

mod twister;
pub use self::twister::Twister;

mod bender;
pub use self::bender::Bender;

mod boolean;
pub use self::boolean::{Intersection, Union};

mod sphere;
pub use self::sphere::Sphere;

mod cylinder;
pub use self::cylinder::{Cone, Cylinder};

mod plane;
pub use self::plane::{NormalPlane, PlaneNegX, PlaneNegY, PlaneNegZ, PlaneX, PlaneY, PlaneZ};

mod mesh;
pub use self::mesh::Mesh;

#[cfg(test)]
mod test;

/// This struct configures evaluation of rounded edges between object.
/// The edge is evaluated in a different more computationally expensive way.
pub struct PrimitiveParameters<S> {
    /// Fade from standard object evaluation to edge evaluation on this fraction of the edge.
    pub fade_range: S,
    /// How much to extend the radius for edge evaluation mode.
    pub r_multiplier: S,
}

const ALWAYS_PRECISE: f32 = 1.;
const EPSILON: f32 = 1e-10;

/// Get a normal from an Object a some point. Do this using approximating the derivative with
/// deltas.
fn normal_from_object<S: Debug + RealField + Float + From<f32>>(
    f: &dyn Object<S>,
    p: &na::Point3<S>,
) -> na::Vector3<S> {
    let null: S = From::from(0.0);
    let e: S = From::from(EPSILON);
    let a: S = From::from(ALWAYS_PRECISE);
    let epsilon_x = na::Vector3::<S>::new(e, null, null);
    let epsilon_y = na::Vector3::<S>::new(null, e, null);
    let epsilon_z = na::Vector3::<S>::new(null, null, e);
    let center = f.approx_value(p, a);
    let dx = f.approx_value(&(p + epsilon_x), a) - center;
    let dy = f.approx_value(&(p + epsilon_y), a) - center;
    let dz = f.approx_value(&(p + epsilon_z), a) - center;
    na::Vector3::<S>::new(dx, dy, dz).normalize()
}

/// Object is the basic trait for any 3d implicit function.
pub trait Object<S: RealField + Float + From<f32>>: ObjectClone<S> + Debug + Sync + Send {
    /// Get the Bounding Box of this Object.
    fn bbox(&self) -> &BoundingBox<S>;
    /// Explicitly set the Bounding Box.
    fn set_bbox(&mut self, _: &BoundingBox<S>) {
        unimplemented!();
    }
    /// Allows to set parameters.
    fn set_parameters(&mut self, _: &PrimitiveParameters<S>) {}
    /// Value is 0 on object surfaces, negative inside and positive outside of objects.
    /// If positive, value is guarateed to be the minimum distance to the object surface.
    /// return some approximation (which is always larger then the proper value).
    /// Only do a proper calculation, for values smaller then slack.
    fn approx_value(&self, _: &na::Point3<S>, _: S) -> S {
        unimplemented!();
    }
    /// Evaluate the normal of ```self``` at the given point.
    fn normal(&self, _: &na::Point3<S>) -> na::Vector3<S> {
        unimplemented!();
    }
    /// Return a translated version of ```self```.
    fn translate(&self, v: &na::Vector3<S>) -> Box<dyn Object<S>> {
        AffineTransformer::new_translate(self.clone_box(), v)
    }
    /// Return a rotated version of ```self```.
    fn rotate(&self, r: &na::Vector3<S>) -> Box<dyn Object<S>> {
        AffineTransformer::new_rotate(self.clone_box(), r)
    }
    /// Return a scaled version of ```self```.
    fn scale(&self, s: &na::Vector3<S>) -> Box<dyn Object<S>> {
        AffineTransformer::new_scale(self.clone_box(), s)
    }
}

/// Trait to allow cloning of ```Box<Object<_>>```.
pub trait ObjectClone<S> {
    /// Clone ```Box<Object<_>>```.
    fn clone_box(&self) -> Box<dyn Object<S>>;
}

impl<S: RealField + Float + From<f32>, T> ObjectClone<S> for T
where
    T: 'static + Object<S> + Clone,
{
    fn clone_box(&self) -> Box<dyn Object<S>> {
        Box::new(self.clone())
    }
}

// We can now implement Clone manually by forwarding to clone_box.
impl<S> Clone for Box<dyn Object<S>> {
    fn clone(&self) -> Box<dyn Object<S>> {
        self.clone_box()
    }
}

// Objects never equal each other
impl<S> PartialEq for Box<dyn Object<S>> {
    fn eq(&self, _: &Box<dyn Object<S>>) -> bool {
        false
    }
}

// Objects are never ordered
impl<S> PartialOrd for Box<dyn Object<S>> {
    fn partial_cmp(&self, _: &Box<dyn Object<S>>) -> Option<::std::cmp::Ordering> {
        None
    }
}
