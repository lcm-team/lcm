# Test overview

The tests are grouped by the behavior they verify.

## Integration tests

| Test | Verifies |
|---|---|
| `test_initial_particle_population` | Initial particle placement and counts |
| `test_particle_exchange_mpi` | Four-rank diagonal exchange and complete state preservation |
| `test_particle_exchange_mpi_edge_cases` | Empty exchange, one-sided traffic, ownership, and conservation |
| `test_lcm_advance_particle_exchange_mpi` | Particle exchange through the full `lcm_advance` path |
| `test_particle_source` | Particle source creation |
| `test_minimal_les_host_*` | Host coupling, cell crossing, and physical boundary behavior |

## Unit tests

| Test | Verifies |
|---|---|
| `test_particle_cell_ownership`, `test_particle_cell_removal` | Local cell ownership and removal |
| `test_particle_migration_metadata` | Metadata updates during migration |
| `test_outgoing_particle_message_packing`, `test_particle_message`, `test_particle_mpi_datatype` | MPI message contents and datatype layout |
| `test_particle_count_exchange`, `test_particle_exchange_buffers` | Exchange counts and buffers |
| `test_exchange_directions`, `test_periodic_particle_destination_wrapping` | Direction selection and periodic destinations |
| `test_particle_location`, `test_local_subdomain_bounds`, `test_vertical_boundary_bounds` | Location and domain-boundary classification |
| `test_velocity_interpolation`, `test_particle_radius` | Velocity interpolation and particle radius |
| `test_parallel_state`, `test_boundary_config`, `test_field_attachment` | Parallel configuration and host-field attachment |

The migration tests cover particle transport only. Two-way thermodynamic coupling will require separate tests.
