module bit_info

use iso_c_binding, only: c_int, c_float, c_double, c_char, c_size_t, c_ptr

implicit none

interface

    subroutine shave_float(a, n_elem, n, shaved, config_file) bind(C, name="shave_float")
      import :: c_float, c_size_t, c_int, c_char
      real(c_float) :: a(*), shaved(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    subroutine halfshave_float(a, n_elem, n, halfshaved, config_file) bind(C, name="halfshave_float")
      import :: c_float, c_size_t, c_int, c_char
      real(c_float) :: a(*), halfshaved(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    subroutine set_float(a, n_elem, n, set, config_file) bind(C, name="set_float")
      import :: c_float, c_size_t, c_int, c_char
      real(c_float) :: a(*), set(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    subroutine groom_float(a, n_elem, n, groomed, config_file) bind(C, name="groom_float")
      import :: c_float, c_size_t, c_int, c_char
      real(c_float) :: a(*), groomed(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    function preserved_information_float(a, b, n_elem, config_file) result(preserved) &
        bind(C, name="preserved_information_float")
      import :: c_float, c_double, c_size_t, c_char
      real(c_float) :: a(*), b(*)
      integer(c_size_t), value :: n_elem
      character(kind=c_char) :: config_file(*)
      real(c_double) :: preserved
    end function

    function pick_bits_to_shave_float(a, n_elem, tolerance, nbits_old, config_file) result(bits) &
        bind(C, name="pick_bits_to_shave_float")
      import :: c_float, c_double, c_size_t, c_int, c_char
      real(c_float) :: a(*)
      integer(c_size_t), value :: n_elem
      real(c_double), value :: tolerance
      integer(c_int), value :: nbits_old
      character(kind=c_char) :: config_file(*)
      integer(c_int) :: bits
    end function

    subroutine shave_double(a, n_elem, n, shaved, config_file) bind(C, name="shave_double")
      import :: c_double, c_size_t, c_int, c_char
      real(c_double) :: a(*), shaved(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    subroutine halfshave_double(a, n_elem, n, halfshaved, config_file) bind(C, name="halfshave_double")
      import :: c_double, c_size_t, c_int, c_char
      real(c_double) :: a(*), halfshaved(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    subroutine bit_rounding_float(a, n_elem, n, rounded, config_file) bind(C, name="bit_rounding_float")
      import :: c_float, c_size_t, c_int, c_char
      real(c_float) :: a(*), rounded(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    subroutine bit_rounding_double(a, n_elem, n, rounded, config_file) bind(C, name="bit_rounding_double")
      import :: c_double, c_size_t, c_int, c_char
      real(c_double) :: a(*), rounded(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    subroutine set_double(a, n_elem, n, set, config_file) bind(C, name="set_double")
      import :: c_double, c_size_t, c_int, c_char
      real(c_double) :: a(*), set(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    subroutine groom_double(a, n_elem, n, groomed, config_file) bind(C, name="groom_double")
      import :: c_double, c_size_t, c_int, c_char
      real(c_double) :: a(*), groomed(*)
      integer(c_size_t), value :: n_elem
      integer(c_int), value :: n
      character(kind=c_char) :: config_file(*)
    end subroutine

    function preserved_information_double(a, b, n_elem, config_file) result(preserved) &
        bind(C, name="preserved_information_double")
      import :: c_double, c_size_t, c_char
      real(c_double) :: a(*), b(*)
      integer(c_size_t), value :: n_elem
      character(kind=c_char) :: config_file(*)
      real(c_double) :: preserved
    end function

    function pick_bits_to_shave_double(a, n_elem, tolerance, nbits_old, config_file) result(bits) &
        bind(C, name="pick_bits_to_shave_double")
      import :: c_double, c_size_t, c_int, c_char
      real(c_double) :: a(*)
      integer(c_size_t), value :: n_elem
      real(c_double), value :: tolerance
      integer(c_int), value :: nbits_old
      character(kind=c_char) :: config_file(*)
      integer(c_int) :: bits
    end function

    function pick_bits_to_shave_binary_search_float(a, n_elem, tolerance, nbits_old, config_file, shaved) result(bits) &
        bind(C, name="pick_bits_to_shave_binary_search_float")
      import :: c_float, c_double, c_size_t, c_int, c_char
      real(c_float) :: a(*), shaved(*)
      integer(c_size_t), value :: n_elem
      real(c_double), value :: tolerance
      integer(c_int), value :: nbits_old
      character(kind=c_char) :: config_file(*)
      integer(c_int) :: bits
    end function

    function pick_bits_to_shave_binary_search_double(a, n_elem, tolerance, nbits_old, config_file, shaved) result(bits) &
        bind(C, name="pick_bits_to_shave_binary_search_double")
      import :: c_double, c_size_t, c_int, c_char
      real(c_double) :: a(*), shaved(*)
      integer(c_size_t), value :: n_elem
      real(c_double), value :: tolerance
      integer(c_int), value :: nbits_old
      character(kind=c_char) :: config_file(*)
      integer(c_int) :: bits
    end function

    function entropy_float(a, n_elem) result(h) bind(C, name="entropy_float")
      import :: c_float, c_double, c_size_t
      real(c_float) :: a(*)
      integer(c_size_t), value :: n_elem
      real(c_double) :: h
    end function

    function entropy_double(a, n_elem) result(h) bind(C, name="entropy_double")
      import :: c_double, c_size_t
      real(c_double) :: a(*)
      integer(c_size_t), value :: n_elem
      real(c_double) :: h
    end function

    function real_info_last_fixed_point_error() result(message) bind(C, name="real_info_last_fixed_point_error")
      import :: c_ptr
      type(c_ptr) :: message
    end function

end interface

interface shave
  procedure :: shave_float, shave_double
end interface

interface halfshave
  procedure :: halfshave_float, halfshave_double
end interface

interface bit_rounding
  procedure :: bit_rounding_float, bit_rounding_double
end interface

interface set
  procedure :: set_float, set_double
end interface

interface groom
  procedure :: groom_float, groom_double
end interface

interface preserved_information
  procedure :: preserved_information_float, preserved_information_double
end interface

interface pick_bits_to_shave
  procedure :: pick_bits_to_shave_float, pick_bits_to_shave_double
end interface

interface pick_bits_to_shave_binary_search
  procedure :: pick_bits_to_shave_binary_search_float, pick_bits_to_shave_binary_search_double
end interface

interface entropy
  procedure :: entropy_float, entropy_double
end interface

end module bit_info
