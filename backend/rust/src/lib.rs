mod climate;
mod distance;
mod generation;
mod eustasy;
mod graph;
mod hydrology;
mod ice_step;
mod lakes;
mod plastic;
mod render;
mod spatial;
mod tectonics;

use pyo3::prelude::*;

#[pymodule]
fn mantle_bloom_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<graph::MarchGraph>()?;
    m.add_class::<spatial::KdIndex>()?;
    m.add_function(wrap_pyfunction!(ice_step::step_ice, m)?)?;
    m.add_function(wrap_pyfunction!(eustasy::seeded_flood_levels, m)?)?;
    m.add_function(wrap_pyfunction!(eustasy::solve_sea_level_connected, m)?)?;
    m.add_function(wrap_pyfunction!(eustasy::below_sea_components, m)?)?;
    m.add_function(wrap_pyfunction!(eustasy::connected_ocean_mask, m)?)?;
    m.add_function(wrap_pyfunction!(climate::humidity_zonal_sweep_rs, m)?)?;
    m.add_function(wrap_pyfunction!(climate::smooth_along_coast_rs, m)?)?;
    m.add_function(wrap_pyfunction!(climate::offset_indices_rs, m)?)?;
    m.add_function(wrap_pyfunction!(climate::moisture_divergence_rs, m)?)?;
    m.add_function(wrap_pyfunction!(climate::smooth_field_rs, m)?)?;
    m.add_function(wrap_pyfunction!(climate::march_substep_prep, m)?)?;
    m.add_function(wrap_pyfunction!(climate::apply_march_abrasion, m)?)?;
    m.add_function(wrap_pyfunction!(render::render_pool_threads, m)?)?;
    m.add_function(wrap_pyfunction!(render::grid_edt, m)?)?;
    m.add_function(wrap_pyfunction!(render::smooth_biome_field_blend_rs, m)?)?;
    m.add_function(wrap_pyfunction!(render::biome_relative_shade_factor_rs, m)?)?;
    m.add_function(wrap_pyfunction!(render::fill_rects_rs, m)?)?;
    m.add_function(wrap_pyfunction!(render::spill_path_rs, m)?)?;
    m.add_function(wrap_pyfunction!(render::render_biome_or_combined, m)?)?;
    m.add_function(wrap_pyfunction!(render::render_elevation, m)?)?;
    m.add_function(wrap_pyfunction!(render::render_constants_rs, m)?)?;
    m.add_function(wrap_pyfunction!(lakes::lake_components_rs, m)?)?;
    m.add_function(wrap_pyfunction!(lakes::catchment_roots_rs, m)?)?;
    m.add_function(wrap_pyfunction!(lakes::lake_hierarchy_events_rs, m)?)?;
    m.add_function(wrap_pyfunction!(lakes::integrate_body_volume_rs, m)?)?;
    m.add_function(wrap_pyfunction!(hydrology::flow_direction_rs, m)?)?;
    m.add_function(wrap_pyfunction!(hydrology::route_downstream_rs, m)?)?;
    m.add_function(wrap_pyfunction!(hydrology::routing_order_rs, m)?)?;
    m.add_function(wrap_pyfunction!(hydrology::sweep_river_reaches_rs, m)?)?;
    m.add_function(wrap_pyfunction!(generation::priority_flood_rs, m)?)?;
    m.add_function(wrap_pyfunction!(generation::nearest_distance_rs, m)?)?;
    m.add_function(wrap_pyfunction!(generation::nearest_seed_labels_rs, m)?)?;
    m.add_function(wrap_pyfunction!(generation::label_component_count_rs, m)?)?;
    m.add_function(wrap_pyfunction!(tectonics::magma_destination_pairs_rs, m)?)?;
    m.add_function(wrap_pyfunction!(tectonics::boundary_nearest_lists_rs, m)?)?;
    m.add_function(wrap_pyfunction!(tectonics::collision_front_labels_rs, m)?)?;
    m.add_function(wrap_pyfunction!(tectonics::gap_contains_rs, m)?)?;
    m.add_function(wrap_pyfunction!(tectonics::magma_allocate_rs, m)?)?;
    Ok(())
}
