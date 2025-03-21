        subroutine setstresstensor(mi, intpt, num_elem, pdata, stx)

            implicit none

            ! io variables
            integer :: mi(*)           ! maximum value of integration point per element (mi[0])
            integer :: num_elem        ! number of elements
            integer :: intpt           ! intercept added to stx
            real(8) :: pdata(*)        ! data received from preCICE
            real(8) :: stx(6, mi(1),*) ! stress vector used by CalculiX

            ! local variables
            integer :: i, j , k, ngp_max, count

            write(*,*) 'setstresstensor: num_elem = ', num_elem
            write(*,*) 'setstresstensor: intpt = ', intpt
            write(*,*) 'setstresstensor: mi(1) = ', mi(1)

            ngp_max = mi(1)

            count = 1
            do i = 1, num_elem
                do j = 1, ngp_max
                    stx(intpt, j, i) = pdata(count)
                    stx(intpt+1, j, i) = pdata(count+1)
                    stx(intpt+2, j, i) = pdata(count+2)
                    count = count + 3
                end do
            end do

        end subroutine setstresstensor
