        SUBROUTINE PRECICE_MULTISCALE_GET_STRAIN()



        END SUBROUTINE PRECICE_MULTISCALE_GET_STRAIN



        SUBROUTINE PRECICE_MULTISCALE_SET_STX(MI,
     &    INTPT, NUM_ELEM, PDATA, STX)

            IMPLICIT NONE

            ! IO VARIABLES
            INTEGER          :: MI(*)                   ! MAXIMUM VALUE OF INTEGRATION POINT PER ELEMENT (MI[0])
            INTEGER          :: NUM_ELEM                ! NUMBER OF ELEMENTS
            INTEGER          :: INTPT                   ! INTERCEPT ADDED TO STX
            REAL(8)          :: PDATA(*)                ! DATA RECEIVED FROM PRECICE
            REAL(8)          :: STX(6, MI(1),*)         ! STRESS VECTOR USED BY CALCULIX


            ! LOCAL VARIABLES
            INTEGER          :: I, J, K, NGP_MAX, COUNT


            NGP_MAX = MI(1)


            COUNT= 1
            DO I = 1, NUM_ELEM
                DO J = 1, NGP_MAX
                    STX(INTPT, J, I) = PDATA(COUNT)
                    STX(INTPT+1, J, I) = PDATA(COUNT+1)
                    STX(INTPT+2, J, I) = PDATA(COUNT+2)
                    COUNT = COUNT + 3
                ENDDO
            ENDDO


        END SUBROUTINE PRECICE_MULTISCALE_SET_STX



        SUBROUTINE PRECICE_MULTISCALE_SET_XSTIFF(MAX_NGP,
     &    INTPT, NUM_ELEM, pdata, XSTIFF)

        IMPLICIT NONE

        ! IO VARIABLES
        INTEGER          :: MAX_NGP       ! MAXIMUM VALUE OF INTEGRATION POINT PER ELEMENT (MI[0])
        INTEGER          :: NUM_ELEM    ! NUMBER OF ELEMENTS
        INTEGER          :: INTPT       ! INTERCEPT ADDED TO XSTIFF
        REAL(8)          :: PDATA(:)    ! DATA RECEIVED FROM PRECICE
        REAL(8)          :: XSTIFF(:)   ! MATERIAL STIFFNESS MATRIX USED BY CALCULIX


        ! LOCAL VARIABLES
        INTEGER          :: I, J, K, NGP_MAX, COUNT


        NGP_MAX = MI(1)

        COUNT= 1
        DO I = 1, NUM_ELEM
            DO J = 1, NGP_MAX
!                XSTIFF(INTPT, J, I) = PDATA(COUNT)
!                XSTIFF(INTPT+1, J, I) = PDATA(COUNT+1)
!                XSTIFF(INTPT+2, J, I) = PDATA(COUNT+2)
                COUNT = COUNT + 3
            ENDDO
        ENDDO

        ! IF(INTPT == 19) THEN
        DO I = 1, NUM_ELEM
            DO J = 1, NGP_MAX
                write(*,*) XSTIFF(1:3, J, I)
            ENDDO
        ENDDO
        ! ENDIF



        END SUBROUTINE PRECICE_MULTISCALE_SET_XSTIFF

