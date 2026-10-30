program test_fortran_abi
  use iso_c_binding, only: c_float, c_double, c_int, c_size_t, c_char, c_null_char, c_associated
  use bit_info
  implicit none

  character(len=1024) :: config_path
  character(kind=c_char) :: config(1025)
  real(c_float) :: single(3), single_out(3)
  real(c_double) :: double(3), double_out(3)
  integer(c_int) :: count
  integer :: index
  integer(c_size_t), parameter :: length = 3_c_size_t

  call get_command_argument(1, config_path)
  if (len_trim(config_path) == 0) error stop 'Missing configuration path'
  config = c_null_char
  do index = 1, len_trim(config_path)
    config(index) = config_path(index:index)
  end do
  single = 2.5_c_float
  double = -2.5_c_double

  single_out = 0
  call shave(single, length, 1_c_int, single_out, config)
  if (any(single_out /= single)) error stop 'shave float'
  double_out = 0
  call shave(double, length, 1_c_int, double_out, config)
  if (any(double_out /= double)) error stop 'shave double'

  single_out = 0
  call halfshave(single, length, 1_c_int, single_out, config)
  if (any(single_out /= single)) error stop 'halfshave float'
  double_out = 0
  call halfshave(double, length, 1_c_int, double_out, config)
  if (any(double_out /= double)) error stop 'halfshave double'

  single_out = 0
  call bit_rounding(single, length, 1_c_int, single_out, config)
  if (any(single_out /= single)) error stop 'round float'
  double_out = 0
  call bit_rounding(double, length, 1_c_int, double_out, config)
  if (any(double_out /= double)) error stop 'round double'

  single_out = 0
  call set(single, length, 1_c_int, single_out, config)
  if (any(single_out /= single)) error stop 'set float'
  double_out = 0
  call set(double, length, 1_c_int, double_out, config)
  if (any(double_out /= double)) error stop 'set double'

  single_out = 0
  call groom(single, length, 1_c_int, single_out, config)
  if (any(single_out /= single)) error stop 'groom float'
  double_out = 0
  call groom(double, length, 1_c_int, double_out, config)
  if (any(double_out /= double)) error stop 'groom double'

  if (preserved_information(single, single_out, length, config) /= 1.0_c_double) &
    error stop 'preserved information float'
  if (preserved_information(double, double_out, length, config) /= 1.0_c_double) &
    error stop 'preserved information double'
  if (pick_bits_to_shave(single, length, 0.99_c_double, 0_c_int, config) /= 0) &
    error stop 'selection float'
  if (pick_bits_to_shave(double, length, 0.99_c_double, 0_c_int, config) /= 0) &
    error stop 'selection double'

  single_out = 0
  count = pick_bits_to_shave_binary_search(single, length, 0.99_c_double, 0_c_int, config, single_out)
  if (count /= 0 .or. any(single_out /= single)) error stop 'binary selection float'
  double_out = 0
  count = pick_bits_to_shave_binary_search(double, length, 0.99_c_double, 0_c_int, config, double_out)
  if (count /= 0 .or. any(double_out /= double)) error stop 'binary selection double'

  if (entropy(single, length) /= 0.0_c_double) error stop 'entropy float'
  if (entropy(double, length) /= 0.0_c_double) error stop 'entropy double'
  if (c_associated(real_info_last_fixed_point_error())) error stop 'unexpected error'

  single = [1.0_c_float, 1.5_c_float, 2.0_c_float]
  double = [1.0_c_double, 1.5_c_double, 2.0_c_double]
  call shave(single, length, 63_c_int, single_out, config)
  if (any(single_out /= single(1))) error stop 'fixed-point range float'
  call shave(double, length, 63_c_int, double_out, config)
  if (any(double_out /= double(1))) error stop 'fixed-point range double'
end program