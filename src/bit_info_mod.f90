module bit_info

use kinds, only: ccs_int
use iso_c_binding, only: c_int, c_float, c_double

implicit none

interface 

  module subroutine confirm(a)
    integer(ccs_int) :: a
  end subroutine

  module subroutine shave_float(a, n_elem, n, shaved) bind(C, name="shave_float")
    real(c_float), dimension(*) :: a
    integer(c_int), value :: n_elem
    integer(c_int), value :: n
    real(c_float), dimension(*) :: shaved
  end subroutine 

  module function preserved_information_float(a, b, n_elem) result(preserved) bind(C, name="preserved_information_float")
    real(c_float), dimension(*) :: a
    real(c_float), dimension(*) :: b
    integer(c_int), value :: n_elem
    real(c_double) :: preserved
  end function 

  module function pick_bits_to_shave_float(a, n_elem, tolerance, nbits_old) result(bits) bind(C, name="pick_bits_to_shave_float")
    real(c_float), dimension(*) :: a
    integer(c_int), value :: n_elem
    real(c_double), value :: tolerance
    integer(c_int), value :: nbits_old
    integer(c_int) :: bits
  end function 

  module subroutine shave_double(a, n_elem, n, shaved) bind(C, name="shave_double")
    real(c_double), dimension(*) :: a
    integer(c_int), value :: n_elem
    integer(c_int), value :: n
    real(c_double), dimension(*) :: shaved
  end subroutine 

  module function preserved_information_double(a, b, n_elem) result(preserved) bind(C, name="preserved_information_double")
    real(c_double), dimension(*) :: a
    real(c_double), dimension(*) :: b
    integer(c_int), value :: n_elem
    real(c_double) :: preserved
  end function 

  module function pick_bits_to_shave_double(a, n_elem, tolerance, nbits_old) result(bits) bind(C, name="pick_bits_to_shave_double")
    real(c_double), dimension(*) :: a
    integer(c_int), value :: n_elem
    real(c_double), value :: tolerance
    integer(c_int), value :: nbits_old
    integer(c_int) :: bits
  end function 

  module function pick_bits_to_shave_binary_search_float(a, n_elem, tolerance, nbits_old) result(bits) bind(C, name="pick_bits_to_shave_binary_search_float")
    real(c_float), dimension(*) :: a
    integer(c_int), value :: n_elem
    real(c_double), value :: tolerance
    integer(c_int), value :: nbits_old
    integer(c_int) :: bits
  end function 

  module function pick_bits_to_shave_binary_search_double(a, n_elem, tolerance, nbits_old) result(bits) bind(C, name="pick_bits_to_shave_binary_search_double")
    real(c_double), dimension(*) :: a
    integer(c_int), value :: n_elem
    real(c_double), value :: tolerance
    integer(c_int), value :: nbits_old
    integer(c_int) :: bits
  end function 

end interface

interface shave
  module procedure shave_float
  module procedure shave_double
end interface

interface preserved_information
  module procedure preserved_information_float
  module procedure preserved_information_double
end interface

interface pick_bits_to_shave
  module procedure pick_bits_to_shave_float
  module procedure pick_bits_to_shave_double
end interface

interface pick_bits_to_shave_binary_search
  module procedure pick_bits_to_shave_binary_search_float
  module procedure pick_bits_to_shave_binary_search_double
end interface

end module bit_info
