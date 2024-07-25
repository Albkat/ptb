program submatrix_method

    use, intrinsic :: iso_fortran_env, only: real64, int64
    implicit none

    real(kind=real64) :: A(4,4), B(4,4)

    A = reshape( (/ 1,0,3,4, &
                    0,4,0,2, &
                    3,0,1,0, &
                    4,2,0,3/), shape(A) )

    print *, "Input matrix:"
    call print_matrix(A)
    print *

    call apply_submatrix_method(A, demo, B)

    print *, "Approximate result matrix:"
    call print_matrix(B)
    print *

    contains

    subroutine apply_submatrix_method(input_matrix, routine, result)
        real(kind=real64), dimension(:,:), intent(in) :: input_matrix
        real(kind=real64), dimension(:,:), intent(out) :: result
        external :: routine
        integer(kind=int64) :: N(2), dim, i, j, k, nonzero, result_col
        real(kind=real64), dimension(:,:), allocatable :: submatrix, submatrix_result
        integer(kind=int64), dimension(:), allocatable :: indices

        interface
            subroutine routine(input_matrix, result)
                use, intrinsic :: iso_fortran_env, only: real64
                real(kind=real64), dimension(:,:), intent(in) :: input_matrix
                real(kind=real64), dimension(:,:), intent(out) :: result
            end subroutine routine
        end interface

        N = shape(input_matrix)
        dim = N(2)

        ! We generate and process a submatrix for each column i of the input matrix
        do i = 1, dim

            ! The size of the submatrix is determined by the number of nonzero elements
            ! in the inducing column i.
            nonzero = count(input_matrix(:,i) /= 0)
            allocate(submatrix(nonzero, nonzero), submatrix_result(nonzero, nonzero), indices(nonzero))


            ! We need to determine which elements in that column are nonzero because
            ! they determine the elements copied into the submatrix.
            k = 1
            do j = 1, dim
                if (input_matrix(j,i) /= 0) then
                    indices(k) = j

                    ! We should take note which column in the submatrix contains elements from the i'th column
                    ! of the original matrix because that will be our result column.
                    ! Note: Elements on the diagonal of the input matrix always have to be nonzero!
                    if (j == i) then
                        result_col = k
                    end if

                    k = k+1
                end if
            end do

            ! Building the submatrix now just means selecting the right rows and columns of the original matrix.
            submatrix = input_matrix(indices,indices)
            print *, "Submatrix", i
            call print_blowed_matrix(nonzero,submatrix)
            print *

            ! Apply the function of interest to the submatrix.
            call routine(submatrix, submatrix_result)
            print *, "Processed submatrix", i
            call print_blowed_matrix(nonzero,submatrix_result)
            print *

            stop

            ! Now we copy the matching result column into the approximate result matrix.
            result(indices,i) = submatrix_result(:,result_col)

            deallocate(submatrix, submatrix_result, indices)
        end do
    end subroutine apply_submatrix_method

    ! A very simple matrix function for demonstration purposes. This could be any matrix function.
    subroutine demo(input_matrix, result)
        real(kind=real64), dimension(:,:), intent(in) :: input_matrix
        real(kind=real64), dimension(:,:), intent(out) :: result
        result = - input_matrix * 2.0d0
    end subroutine demo

    subroutine print_matrix(input_matrix)
        real(kind=real64), dimension(:,:), intent(in) :: input_matrix
        integer(kind=int64) :: N(2), i
        N = shape(input_matrix, kind=int64)
        do i = 1, N(1)
            print *, input_matrix(i,:)
        end do
    end subroutine print_matrix

end program submatrix_method