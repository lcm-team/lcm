!> Flat directional message buffers for standalone particle MPI exchange.
module lcm_particle_exchange_buffers
   use lcm_exchange_directions, only: lcm_exchange_direction_count, &
                                      lcm_exchange_direction_for_phase, &
                                      lcm_exchange_directions_per_phase, &
                                      lcm_require_valid_exchange_phase
   use lcm_types, only: lcm_particle_mpi_message_t

   implicit none

   private

   public :: lcm_clear_particle_exchange_buffers
   public :: lcm_clear_particle_exchange_phase_buffers
   public :: lcm_append_send_particle_message
   public :: lcm_append_receive_particle_message
   public :: lcm_send_particle_message_count
   public :: lcm_receive_particle_message_count
   public :: lcm_expected_receive_particle_message_count
   public :: lcm_reserve_receive_particle_messages
   public :: lcm_get_send_particle_message_view
   public :: lcm_get_receive_particle_message_view
   public :: lcm_mark_receive_particle_messages_received
   public :: lcm_get_send_particle_message
   public :: lcm_get_receive_particle_message

   type :: lcm_particle_message_buffer_t
      type(lcm_particle_mpi_message_t), allocatable :: messages(:)
      integer :: count = 0
      integer :: expected_count = 0
   end type lcm_particle_message_buffer_t

   type(lcm_particle_message_buffer_t), target, save :: &
      send_buffers(lcm_exchange_direction_count)
   type(lcm_particle_message_buffer_t), target, save :: &
      receive_buffers(lcm_exchange_direction_count)

contains

   !> Reset message counts while retaining allocated storage for the next step.
   subroutine lcm_clear_particle_exchange_buffers()
      integer :: direction

      do direction = 1, lcm_exchange_direction_count
         send_buffers(direction)%count = 0
         receive_buffers(direction)%count = 0
         receive_buffers(direction)%expected_count = 0
      enddo
   end subroutine lcm_clear_particle_exchange_buffers

   !> Reset only the two buffers used by one cardinal exchange phase.
   subroutine lcm_clear_particle_exchange_phase_buffers(phase)
      integer, intent(in) :: phase
      integer :: direction
      integer :: phase_direction

      call lcm_require_valid_exchange_phase(phase)
      do phase_direction = 1, lcm_exchange_directions_per_phase
         direction = lcm_exchange_direction_for_phase(phase, phase_direction)
         send_buffers(direction)%count = 0
         receive_buffers(direction)%count = 0
         receive_buffers(direction)%expected_count = 0
      enddo
   end subroutine lcm_clear_particle_exchange_phase_buffers

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

   !> Number of messages reported by the matching neighbor during count exchange.
   !> This reserves receive storage but does not mark payload messages as received.
   integer function lcm_expected_receive_particle_message_count(direction) result(count)
      integer, intent(in) :: direction

      call require_valid_direction(direction)
      count = receive_buffers(direction)%expected_count
   end function lcm_expected_receive_particle_message_count

   subroutine lcm_reserve_receive_particle_messages(direction, expected_count)
      integer, intent(in) :: direction
      integer, intent(in) :: expected_count

      call require_valid_direction(direction)
      if (expected_count < 0) error stop 'LCM expected receive count must not be negative'

      receive_buffers(direction)%expected_count = expected_count
      call ensure_buffer_capacity(receive_buffers(direction), expected_count)
   end subroutine lcm_reserve_receive_particle_messages

   !> Provide direct contiguous access to the active send-buffer payload.
   !>
   !> The view is valid only until the matching direction is cleared or grown.
   !> Transport must not modify the messages or retain the view after it waits.
   subroutine lcm_get_send_particle_message_view(direction, messages)
      integer, intent(in) :: direction
      type(lcm_particle_mpi_message_t), contiguous, pointer, intent(out) :: messages(:)

      call require_valid_direction(direction)
      if (send_buffers(direction)%count == 0) then
         nullify(messages)
      else
         messages => send_buffers(direction)%messages( &
            1:send_buffers(direction)%count)
      endif
   end subroutine lcm_get_send_particle_message_view

   !> Provide direct contiguous access to reserved receive-buffer storage.
   !>
   !> The returned region has the count established by the preceding count
   !> exchange. Transport fills it in place and then marks it as received.
   subroutine lcm_get_receive_particle_message_view(direction, messages)
      integer, intent(in) :: direction
      type(lcm_particle_mpi_message_t), contiguous, pointer, intent(out) :: messages(:)

      call require_valid_direction(direction)
      if (receive_buffers(direction)%expected_count == 0) then
         nullify(messages)
      else
         messages => receive_buffers(direction)%messages( &
            1:receive_buffers(direction)%expected_count)
      endif
   end subroutine lcm_get_receive_particle_message_view

   subroutine lcm_mark_receive_particle_messages_received(direction)
      integer, intent(in) :: direction

      call require_valid_direction(direction)
      receive_buffers(direction)%count = receive_buffers(direction)%expected_count
   end subroutine lcm_mark_receive_particle_messages_received

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

      if (required_capacity <= 0) return
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
