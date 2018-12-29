use alga::general::Real;
use alga::linear::Transformation;
use na;
use num_traits::Float;
use {BoundingBox, Object, PrimitiveParameters};

#[derive(Clone, Debug)]
/// AffineTransformer is a primitive that takes an object as input and allows to modify it using
/// affine transforms.
/// Usually it is used indirectly through ```Object::scale()```, ```Object::translate()``` or ```Object::rotate()```.
pub struct AffineTransformer<S: Real> {
    object: Box<Object<S>>,
    transform: na::Matrix4<S>,
    transposed3x3: na::Matrix3<S>,
    scale_min: S,
    bbox: BoundingBox<S>,
}

impl<S: Real + Float + From<f32>> Object<S> for AffineTransformer<S> {
    fn approx_value(&self, p: &na::Point3<S>, slack: S) -> S {
        let approx = self.bbox.distance(p);
        if approx <= slack {
            self.object
                .approx_value(&self.transform.transform_point(&p), slack / self.scale_min)
                * self.scale_min
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
        let normal_at_p = self.object.normal(&self.transform.transform_point(&p));
        let transformed_normal = self.transposed3x3 * normal_at_p;
        transformed_normal.normalize()
    }
    fn translate(&self, v: &na::Vector3<S>) -> Box<Object<S>> {
        let new_trans = self.transform.append_translation(&-v);
        AffineTransformer::new_with_scaler(self.object.clone(), new_trans, self.scale_min)
    }
    fn rotate(&self, r: &na::Vector3<S>) -> Box<Object<S>> {
        let euler = ::na::Rotation::from_euler_angles(r.x, r.y, r.z).to_homogeneous();
        let new_trans = self.transform * euler;
        AffineTransformer::new_with_scaler(self.object.clone(), new_trans, self.scale_min)
    }
    fn scale(&self, s: &na::Vector3<S>) -> Box<Object<S>> {
        let one: S = From::from(1f32);
        let new_trans = self.transform.append_nonuniform_scaling(&na::Vector3::new(
            one / s.x,
            one / s.y,
            one / s.z,
        ));
        AffineTransformer::new_with_scaler(
            self.object.clone(),
            new_trans,
            self.scale_min * Float::min(s.x, Float::min(s.y, s.z)),
        )
    }
}

impl<S: Real + Float + From<f32>> AffineTransformer<S> {
    fn identity(o: Box<Object<S>>) -> Box<Object<S>> {
        AffineTransformer::new(o, na::Matrix4::identity())
    }
    fn new(o: Box<Object<S>>, t: na::Matrix4<S>) -> Box<AffineTransformer<S>> {
        let one: S = From::from(1f32);
        AffineTransformer::new_with_scaler(o, t, one)
    }
    fn new_with_scaler(
        o: Box<Object<S>>,
        t: na::Matrix4<S>,
        scale_min: S,
    ) -> Box<AffineTransformer<S>> {
        // TODO: Calculate scale_min from t.
        // This should be something similar to
        // 1./Vector::new(t.x.x, t.y.x, t.z.x).magnitude().min(
        // 1./Vector::new(t.x.y, t.y.y, t.z.y).magnitude().min(
        // 1./Vector::new(t.x.z, t.y.z, t.z.z).magnitude()))

        match t.try_inverse() {
            None => panic!("Failed to invert {:?}", t),
            Some(t_inv) => {
                let bbox = o.bbox().transform(&t_inv);
                let transposed3x3 = t
                    .fixed_slice::<::na::core::dimension::U3, ::na::core::dimension::U3>(0, 0)
                    .transpose();
                Box::new(AffineTransformer {
                    object: o,
                    transform: t,
                    transposed3x3,
                    scale_min,
                    bbox,
                })
            }
        }
    }
    /// Create a new translated version of the input.
    pub fn new_translate(o: Box<Object<S>>, v: &na::Vector3<S>) -> Box<Object<S>> {
        AffineTransformer::identity(o).translate(v)
    }
    /// Create a new rotated version of the input.
    pub fn new_rotate(o: Box<Object<S>>, r: &na::Vector3<S>) -> Box<Object<S>> {
        AffineTransformer::identity(o).rotate(r)
    }
    /// Create a new scaled version of the input.
    pub fn new_scale(o: Box<Object<S>>, s: &na::Vector3<S>) -> Box<Object<S>> {
        AffineTransformer::identity(o).scale(s)
    }
}

#[cfg(test)]
mod test {
    use super::super::test::MockObject;
    use super::*;
    use std::sync::mpsc::sync_channel;

    #[test]
    fn translate() {
        let (sender, receiver) = sync_channel(1);
        let normal = na::Vector3::new(1.0, 0.0, 0.0);
        let mock_object = MockObject::new(1.0, normal, sender);
        let translation = na::Vector3::new(0.0001, 0.0, 0.0);
        let translated = mock_object.translate(&translation);
        let p = na::Point3::new(1.0, 0.0, 0.0);
        assert_eq!(translated.normal(&p), normal);
        assert_eq!(receiver.recv().unwrap(), p - translation);
    }

    #[test]
    fn rotate() {
        let (sender, receiver) = sync_channel(1);
        let normal = na::Vector3::new(1.0, 0.0, 0.0);
        let mock_object = MockObject::new(1.0, normal, sender);
        let rotation = na::Vector3::new(0.0, 0.0, ::std::f64::consts::PI / 6.0);
        let rotated = mock_object.rotate(&rotation);
        let p = na::Point3::new(1.0, 0.0, 0.0);

        assert_relative_eq!(
            rotated.normal(&p),
            na::Vector3::new(num_traits::Float::sqrt(3.0) / 2.0, -0.5, 0.0)
        );
        assert_relative_eq!(
            receiver.try_recv().unwrap(),
            na::Point3::new(num_traits::Float::sqrt(3.0) / 2.0, 0.5, 0.0)
        );
    }
}
