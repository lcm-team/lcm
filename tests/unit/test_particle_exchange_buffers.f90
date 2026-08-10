program test_particle_exchange_buffers
   use iso_fortran_env, only: int64
   use lcm_exchange_directions, only: lcm_exchange_east, lcm_exchange_north, &
                                      lcm_exchange_northeast
   use lcm_particle_exchange_buffers, only: lcm_append_receive_particle_message, &
                                            lcm_append_send_particle_message, &
                                            lcm_clear_particle_exchange_buffers, &
                                            lcm_get_receive_particle_message, &
                                            lcm_get_send_particle_message, &
                                            lcm_receive_particle_message_count, &
                                            lcm_send_particle_message_count
   use lcm_types, only: lcm_particle_mpi_message_t

   implicit none

   type(lcm_particle_mpi_message_t) :: message
   type(lcm_particle_mpi_message_t) :: copied_message

   call lcm_clear_particle_exchange_buffers()
   call populate_message(message, 101_int64, 11, 7, 2)
   call lcm_append_send_particle_message(lcm_exchange_east, message)
   call populate_message(message, 102_int64, 12, 7, 2)
   call lcm_append_send_particle_message(lcm_exchange_east, message)
   call populate_message(message, 103_int64, 11, 8, 2)
   call lcm_append_receive_particle_message(lcm_exchange_north, message)

   call assert_equal(lcm_send_particle_message_count(lcm_exchange_east), 2, &
                     'east send count')
   call assert_equal(lcm_send_particle_message_count(lcm_exchange_northeast), 0, &
                     'empty northeast send count')
   call assert_equal(lcm_receive_particle_message_count(lcm_exchange_north), 1, &
                     'north receive count')

   call lcm_get_send_particle_message(lcm_exchange_east, 2, copied_message)
   call assert_equal(int(copied_message%particle%id), 102, 'second east id')
   call assert_equal(copied_message%global_i, 12, 'second east global i')
   call lcm_get_receive_particle_message(lcm_exchange_north, 1, copied_message)
   call assert_equal(int(copied_message%particle%id), 103, 'north receive id')

   call lcm_clear_particle_exchange_buffers()
   call assert_equal(lcm_send_particle_message_count(lcm_exchange_east), 0, &
                     'cleared east send count')
   call assert_equal(lcm_receive_particle_message_count(lcm_exchange_north), 0, &
                     'cleared north receive count')

contains

   subroutine populate_message(message, id, global_i, global_j, k)
      type(lcm_particle_mpi_message_t), intent(out) :: message
      integer(int64), intent(in) :: id
      integer, intent(in) :: global_i
      integer, intent(in) :: global_j
      integer, intent(in) :: k

      message = lcm_particle_mpi_message_t()
      message%global_i = global_i
      message%global_j = global_j
      message%k = k
      message%particle%id = id
   end subroutine populate_message

   subroutine assert_equal(actual, expected, description)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

end program test_particle_exchange_buffers
