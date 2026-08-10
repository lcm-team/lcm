!> Flat directional message buffers for standalone particle MPI exchange.
module lcm_particle_exchange_buffers
   use lcm_exchange_directions, only: lcm_exchange_direction_count
   use lcm_types, only: lcm_particle_mpi_message_t

   implicit none

   private

   public :: lcm_clear_particle_exchange_buffers
   public :: lcm_append_send_particle_message
   public :: lcm_append_receive_particle_message
   public :: lcm_send_particle_message_count
   public :: lcm_receive_particle_message_count
   public :: lcm_get_send_particle_message
   public :: lcm_get_receive_particle_message

   type :: lcm_particle_message_buffer_t
      type(lcm_particle_mpi_message_t), allocatable :: messages(:)
      integer :: count = 0
   end type lcm_particle_message_buffer_t

   type(lcm_particle_message_buffer_t), save :: send_buffers(lcm_exchange_direction_count)
   type(lcm_particle_message_buffer_t), save :: receive_buffers(lcm_exchange_direction_count)

contains

   !> Reset message counts while retaining allocated storage for the next step.
   subroutine lcm_clear_particle_exchange_buffers()
      integer :: direction

      do direction = 1, lcm_exchange_direction_count
         send_buffers(direction)%count = 0
         receive_buffers(direction)%count = 0
      enddo
   end subroutine lcm_clear_particle_exchange_buffers

   subroutine lcm_append_send_particle_message(direction, message)
      integer, intent(in) :: direction
      type(lcm_particle_mpi_message_t), intent(in) :: message

      call append_particle_message(send_buffers, direction, message)
   end subroutine lcm_append_send_particle_message

   subroutine lcm_append_receive_particle_message(direction, message)
      integer, intent(in) :: direction
      type(lcm_particle_mpi_message_t), intent(in) :: message

      call append_particle_message(receive_buffers, direction, message)
   end subroutine lcm_append_receive_particle_message

   integer function lcm_send_particle_message_count(direction) result(count)
      integer, intent(in) :: direction

      call require_valid_direction(direction)
      count = send_buffers(direction)%count
   end function lcm_send_particle_message_count

   integer function lcm_receive_particle_message_count(direction) result(count)
      integer, intent(in) :: direction

      call require_valid_direction(direction)
      count = receive_buffers(direction)%count
   end function lcm_receive_particle_message_count

   subroutine lcm_get_send_particle_message(direction, index, message)
      integer, intent(in) :: direction
      integer, intent(in) :: index
      type(lcm_particle_mpi_message_t), intent(out) :: message

      call get_particle_message(send_buffers, direction, index, message)
   end subroutine lcm_get_send_particle_message

   subroutine lcm_get_receive_particle_message(direction, index, message)
      integer, intent(in) :: direction
      integer, intent(in) :: index
      type(lcm_particle_mpi_message_t), intent(out) :: message

      call get_particle_message(receive_buffers, direction, index, message)
   end subroutine lcm_get_receive_particle_message

   subroutine append_particle_message(buffers, direction, message)
      type(lcm_particle_message_buffer_t), intent(inout) :: buffers(:)
      integer, intent(in) :: direction
      type(lcm_particle_mpi_message_t), intent(in) :: message

      call require_valid_direction(direction)
      call ensure_buffer_capacity(buffers(direction), buffers(direction)%count + 1)

      buffers(direction)%count = buffers(direction)%count + 1
      buffers(direction)%messages(buffers(direction)%count) = message
   end subroutine append_particle_message

   subroutine get_particle_message(buffers, direction, index, message)
      type(lcm_particle_message_buffer_t), intent(in) :: buffers(:)
      integer, intent(in) :: direction
      integer, intent(in) :: index
      type(lcm_particle_mpi_message_t), intent(out) :: message

      call require_valid_direction(direction)
      if (index < 1 .or. index > buffers(direction)%count) then
         error stop 'LCM particle exchange message index is out of range'
      endif

      message = buffers(direction)%messages(index)
   end subroutine get_particle_message

   subroutine ensure_buffer_capacity(buffer, required_capacity)
      type(lcm_particle_message_buffer_t), intent(inout) :: buffer
      integer, intent(in) :: required_capacity
      type(lcm_particle_mpi_message_t), allocatable :: previous_messages(:)
      integer :: new_capacity

      if (.not. allocated(buffer%messages)) then
         allocate(buffer%messages(max(1, required_capacity)))
         return
      endif
      if (required_capacity <= size(buffer%messages)) return

      new_capacity = max(required_capacity, 2 * size(buffer%messages))
      call move_alloc(buffer%messages, previous_messages)
      allocate(buffer%messages(new_capacity))
      buffer%messages(1:buffer%count) = previous_messages(1:buffer%count)
   end subroutine ensure_buffer_capacity

   subroutine require_valid_direction(direction)
      integer, intent(in) :: direction

      if (direction < 1 .or. direction > lcm_exchange_direction_count) then
         error stop 'LCM particle exchange direction is invalid'
      endif
   end subroutine require_valid_direction

end module lcm_particle_exchange_buffers
