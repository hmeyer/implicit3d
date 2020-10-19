use crate::{Object, RealField};
use bbox::BoundingBox;
use num_traits::Float;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};

#[derive(Clone, Debug)]
pub struct MockObject<S: RealField> {
    value: S,
    normal: na::Vector3<S>,
    bbox: BoundingBox<S>,
    normal_call_sender: Option<SyncSender<na::Point3<S>>>,
}

impl<S: ::std::fmt::Debug + Float + RealField> MockObject<S> {
    pub fn new(value: S, normal: na::Vector3<S>) -> Self {
        Self::new_with_bbox(value, normal, BoundingBox::infinity())
    }
    pub fn new_with_bbox(value: S, normal: na::Vector3<S>, bbox: BoundingBox<S>) -> Self {
        MockObject {
            value,
            normal,
            bbox,
            normal_call_sender: None,
        }
    }
    pub fn add_normal_call_recorder(&mut self, buffer_size: usize) -> Receiver<na::Point3<S>> {
        let (sender, receiver) = sync_channel(buffer_size);
        self.normal_call_sender = Some(sender);
        receiver
    }
}

impl<S: ::std::fmt::Debug + RealField + Float + From<f32>> Object<S> for MockObject<S> {
    fn approx_value(&self, _: &na::Point3<S>, _: S) -> S {
        self.value
    }
    fn normal(&self, p: &na::Point3<S>) -> na::Vector3<S> {
        if let Some(sender) = &self.normal_call_sender {
            sender.send(*p).unwrap();
        }
        self.normal
    }
    fn bbox(&self) -> &BoundingBox<S> {
        &self.bbox
    }
}
