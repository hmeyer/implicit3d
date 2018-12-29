use super::Object;
use alga::general::Real;
use bbox::BoundingBox;
use na;
use num_traits::Float;
use std::sync::mpsc::SyncSender;

#[derive(Clone, Debug)]
pub struct MockObject<S: Real> {
    value: S,
    normal: na::Vector3<S>,
    bbox: BoundingBox<S>,
    normal_call_sender: SyncSender<na::Point3<S>>,
}

impl<S: ::std::fmt::Debug + Float + Real> MockObject<S> {
    pub fn new(
        value: S,
        normal: na::Vector3<S>,
        normal_call_sender: SyncSender<na::Point3<S>>,
    ) -> Box<MockObject<S>> {
        Box::new(MockObject {
            value,
            normal,
            bbox: BoundingBox::infinity(),
            normal_call_sender,
        })
    }
}

impl<S: ::std::fmt::Debug + Real + Float + From<f32>> Object<S> for MockObject<S> {
    fn approx_value(&self, _: &na::Point3<S>, _: S) -> S {
        self.value
    }
    fn normal(&self, p: &na::Point3<S>) -> na::Vector3<S> {
        self.normal_call_sender.send(*p).unwrap();
        self.normal
    }
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
}
