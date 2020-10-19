#[macro_use]
extern crate bencher;

use bencher::Bencher;
use implicit3d::{
    Intersection, Object, PlaneNegX, PlaneNegY, PlaneNegZ, PlaneX, PlaneY, PlaneZ, RealField,
    Sphere, Twister,
};
use nalgebra as na;
use num_traits::{Float, FloatConst};
use std::fmt::Debug;

const STEPS: usize = 50;

fn evaluate<S: From<f32> + Debug + Float + RealField>(obj: &dyn Object<S>) -> S {
    let zero = From::from(0f32);
    let mut p = na::Point3::new(zero, zero, obj.bbox().min.z);
    let xd = (obj.bbox().max.x - obj.bbox().min.x) / From::from(STEPS as f32);
    let yd = (obj.bbox().max.y - obj.bbox().min.y) / From::from(STEPS as f32);
    let zd = (obj.bbox().max.z - obj.bbox().min.z) / From::from(STEPS as f32);
    let slack = Float::min(xd, Float::min(yd, zd)) / From::from(10f32);
    let mut result = zero;
    for _ in 0..STEPS {
        p.y = obj.bbox().min.y;
        for _ in 0..STEPS {
            p.x = obj.bbox().min.x;
            for _ in 0..STEPS {
                result += obj.approx_value(&p, slack);
                p.x += xd;
            }
            p.y += yd;
        }
        p.z += zd;
    }
    result
}

fn normals<S: 'static + From<f32> + Debug + Float + RealField>(
    obj: &dyn Object<S>,
) -> na::Vector3<S> {
    let zero = From::from(0f32);
    let mut p = na::Point3::new(zero, zero, obj.bbox().min.z);
    let xd = (obj.bbox().max.x - obj.bbox().min.x) / From::from(STEPS as f32);
    let yd = (obj.bbox().max.y - obj.bbox().min.y) / From::from(STEPS as f32);
    let zd = (obj.bbox().max.z - obj.bbox().min.z) / From::from(STEPS as f32);
    let mut result = na::Vector3::new(zero, zero, zero);
    for _ in 0..STEPS {
        p.y = obj.bbox().min.y;
        for _ in 0..STEPS {
            p.x = obj.bbox().min.x;
            for _ in 0..STEPS {
                result += obj.normal(&p);
                p.x += xd;
            }
            p.y += yd;
        }
        p.z += zd;
    }
    result
}

fn sphere<S: From<f32> + Debug + Float + RealField>(b: &mut Bencher) {
    let object = Sphere::new(From::from(1f32));
    b.iter(|| evaluate(&object as &dyn Object<S>));
}
fn sphere_normals<S: From<f32> + Debug + Float + RealField>(b: &mut Bencher) {
    let object = Sphere::new(From::from(1f32));
    b.iter(|| normals(&object as &dyn Object<S>));
}

fn create_cube<S: From<f32> + Debug + Float + RealField>() -> Box<dyn Object<S>> {
    let zero = From::from(0f32);
    let point_five = From::from(0.5f32);
    Intersection::from_vec(
        vec![
            Box::new(PlaneX::new(point_five)),
            Box::new(PlaneNegX::new(point_five)),
            Box::new(PlaneY::new(point_five)),
            Box::new(PlaneNegY::new(point_five)),
            Box::new(PlaneZ::new(point_five)),
            Box::new(PlaneNegZ::new(point_five)),
        ],
        zero,
    )
    .unwrap()
}

fn cube<S: From<f32> + Debug + Float + RealField>(b: &mut Bencher) {
    let object = create_cube();
    b.iter(|| evaluate(&*object as &dyn Object<S>));
}
fn cube_normals<S: From<f32> + Debug + Float + RealField>(b: &mut Bencher) {
    let object = create_cube();
    b.iter(|| normals(&*object as &dyn Object<S>));
}

fn create_hollow_cube<S: From<f32> + Debug + Float + FloatConst + RealField>() -> Box<dyn Object<S>>
{
    Intersection::difference_from_vec(
        vec![create_cube(), Box::new(Sphere::new(From::from(0.5f32)))],
        From::from(0.2f32),
    )
    .unwrap()
}

fn hollow_cube<S: From<f32> + Debug + Float + FloatConst + RealField>(b: &mut Bencher) {
    let object = create_hollow_cube();
    b.iter(|| evaluate(&*object as &dyn Object<S>));
}
fn hollow_cube_normals<S: From<f32> + Debug + Float + FloatConst + RealField>(b: &mut Bencher) {
    let object = create_hollow_cube();
    b.iter(|| normals(&*object as &dyn Object<S>));
}

fn twisted_cube<S: From<f32> + Debug + Float + FloatConst + RealField>(b: &mut Bencher) {
    let object = Twister::new(create_cube(), From::from(4f32));
    b.iter(|| evaluate(&object as &dyn Object<S>));
}
fn twisted_cube_normals<S: From<f32> + Debug + Float + FloatConst + RealField>(b: &mut Bencher) {
    let object = Twister::new(create_cube(), From::from(4f32));
    b.iter(|| normals(&object as &dyn Object<S>));
}

benchmark_group!(
    bench_values_f32,
    sphere<f32>,
    cube<f32>,
    hollow_cube<f32>,
    twisted_cube<f32>
);
benchmark_group!(
    bench_values_f64,
    sphere<f64>,
    cube<f64>,
    hollow_cube<f64>,
    twisted_cube<f64>
);
benchmark_group!(
    bench_normals_f32,
    sphere_normals<f32>,
    cube_normals<f32>,
    hollow_cube_normals<f32>,
    twisted_cube_normals<f32>
);
benchmark_group!(
    bench_normals_f64,
    sphere_normals<f64>,
    cube_normals<f64>,
    hollow_cube_normals<f64>,
    twisted_cube_normals<f64>
);
benchmark_main!(
    bench_values_f32,
    bench_normals_f32,
    bench_values_f64,
    bench_normals_f64
);
