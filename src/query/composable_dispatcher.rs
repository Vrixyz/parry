//! Logic to dispatch queries to the appropriate handler.
//!
//! This can be customized to add support to your own types.

use crate::{
    math::{Isometry, Real},
    query::QueryDispatcher,
    shape::{Shape, ShapeType},
};

use super::{ContactManifold, PersistentQueryDispatcher};

pub struct DispatchDefinition<ManifoldData, ContactData> {
    pub type_arg1: ShapeType,
    pub type_arg2: ShapeType,
    pub func: Box<
        dyn Fn(
            &Isometry<Real>,
            &dyn Shape,
            &dyn Shape,
            Real,
            &mut ContactManifold<ManifoldData, ContactData>,
        ),
    >,
}

pub struct ComposableQueryDispatcher {
    dispatchers: Vec<((ShapeType, ShapeType), Vec<Box<dyn QueryDispatcher>>)>,
}

impl QueryDispatcher for ComposableQueryDispatcher {
    fn intersection_test(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<bool, super::Unsupported> {
        todo!()
    }

    fn distance(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<Real, super::Unsupported> {
        todo!()
    }

    fn contact(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
    ) -> Result<Option<super::Contact>, super::Unsupported> {
        todo!()
    }

    fn closest_points(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_dist: Real,
    ) -> Result<super::ClosestPoints, super::Unsupported> {
        todo!()
    }

    fn cast_shapes(
        &self,
        pos12: &Isometry<Real>,
        local_vel12: &crate::math::Vector<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        options: super::ShapeCastOptions,
    ) -> Result<Option<super::ShapeCastHit>, super::Unsupported> {
        todo!()
    }

    fn cast_shapes_nonlinear(
        &self,
        motion1: &super::NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &super::NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Result<Option<super::ShapeCastHit>, super::Unsupported> {
        todo!()
    }
}
