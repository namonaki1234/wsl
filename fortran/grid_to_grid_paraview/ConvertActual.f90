!***********************************************************************************************************************************
!***********************************************************************************************************************************
!**** UPACS�̗L�������c�[��                                                                                                     ****
!****                                                                                                                           ****
!**** Name          : ConvertActual.f90                                                                                         ****
!**** Date          : 2019.10.03                                                                                                ****
!**** Programed by  : MANABU UENO                                                                                               ****
!**** Discription   :                                                                                                           ****
!***********************************************************************************************************************************
!***********************************************************************************************************************************
program Main
! �ϐ��錾 *************************************************************************************************************************
implicit none

integer i,j,k,m,n
integer nblocks !�}���`�u���b�N��
double precision mach,alpha,re,t !Phys.q�̃w�b�_�[

double precision,parameter :: rho0 = 7.718858207e+02 !###�v�m�F### input_impg.txt�ɍ��킹�ĕύX cRef
double precision,parameter :: c0 = 3.723909451e+00 !###�v�m�F### input_impg.txt�ɍ��킹�ĕύX roRef
double precision,parameter :: L0 = 1.0e-3 !###�v�m�F### input_impg.txt�ɍ��킹�ĕύX alenRef

type :: PhysData
    double precision,allocatable :: q(:,:,:,:) !�����ʃf�[�^�@������i,j,k,������
    integer :: max(3) !1��imax, 2��jmax, 3��kmax�̒l
end type

type(PhysData),allocatable :: Blk_phys(:)

type :: GridData
    double precision,allocatable :: x(:,:,:)
    double precision,allocatable :: y(:,:,:)
    double precision,allocatable :: z(:,:,:)
    integer :: max(3) !1��imax, 2��jmax, 3��kmax�̒l
end type

type(GridData),allocatable :: Blk_grid(:)

! �����J�n *************************************************************************************************************************
write (*,*) 'Start ConvertGrid'
call ConvertGrid
! write (*,*) 'Start ConvertPhys'
! call ConvertPhys
write (*,*) 'Success!'
! �����I�� *************************************************************************************************************************

! �����葱�� ***********************************************************************************************************************
contains
!***********************************************************************************************************************************
!**** grid.g�̗L������                                                                                                          ****
!***********************************************************************************************************************************
subroutine ConvertGrid

! �����J�n *************************************************************************************************************************
! �t�@�C������**********************************************************************************************************************
    open (unit=7,form='unformatted',file='grid.g')
    read (7) nblocks
    !�u���b�N���̏o��
    print *, nblocks
    allocate (Blk_grid(nblocks))
    read (7) (Blk_grid(m)%max(1),Blk_grid(m)%max(2),Blk_grid(m)%max(3),m=1,nblocks)

    do m = 1,nblocks
        allocate (Blk_grid(m)%x(Blk_grid(m)%max(1),Blk_grid(m)%max(2),Blk_grid(m)%max(3)))
        allocate (Blk_grid(m)%y(Blk_grid(m)%max(1),Blk_grid(m)%max(2),Blk_grid(m)%max(3)))
        allocate (Blk_grid(m)%z(Blk_grid(m)%max(1),Blk_grid(m)%max(2),Blk_grid(m)%max(3)))
        read (7) &
            (((Blk_grid(m)%x(i,j,k),i=1,Blk_grid(m)%max(1)),j=1,Blk_grid(m)%max(2)),k=1,Blk_grid(m)%max(3)), &
            (((Blk_grid(m)%y(i,j,k),i=1,Blk_grid(m)%max(1)),j=1,Blk_grid(m)%max(2)),k=1,Blk_grid(m)%max(3)), &
            (((Blk_grid(m)%z(i,j,k),i=1,Blk_grid(m)%max(1)),j=1,Blk_grid(m)%max(2)),k=1,Blk_grid(m)%max(3))
    end do
    close (7)

! �L������**********************************************************************************************************************
    do m = 1,nblocks
        do k = 1,Blk_grid(m)%max(3)
            do j = 1,Blk_grid(m)%max(2)
                do i = 1,Blk_grid(m)%max(1)
                    Blk_grid(m)%x(i,j,k) = Blk_grid(m)%x(i,j,k)*L0
                    Blk_grid(m)%y(i,j,k) = Blk_grid(m)%y(i,j,k)*L0
                    Blk_grid(m)%z(i,j,k) = Blk_grid(m)%z(i,j,k)*L0
                end do
            end do
        end do
    end do

! �t�@�C���o��**********************************************************************************************************************
    open (unit=7,form='unformatted',file='grid_Paraview.g')
    write (7) nblocks
    write (7) (Blk_grid(m)%max(1),Blk_grid(m)%max(2),Blk_grid(m)%max(3),m=1,nblocks)
    do m = 1,nblocks
        write (7) &
            (((Blk_grid(m)%x(i,j,k),i=1,Blk_grid(m)%max(1)),j=1,Blk_grid(m)%max(2)),k=1,Blk_grid(m)%max(3)), &
            (((Blk_grid(m)%y(i,j,k),i=1,Blk_grid(m)%max(1)),j=1,Blk_grid(m)%max(2)),k=1,Blk_grid(m)%max(3)), &
            (((Blk_grid(m)%z(i,j,k),i=1,Blk_grid(m)%max(1)),j=1,Blk_grid(m)%max(2)),k=1,Blk_grid(m)%max(3))
    end do
    close (7)

    deallocate (Blk_grid)

! �����I�� *************************************************************************************************************************
end subroutine ConvertGrid

!************************************************************************************************************************************
!**** Phys.q�̗L������                                                                                                          ****
!***********************************************************************************************************************************
subroutine ConvertPhys

! �����J�n *************************************************************************************************************************
! �t�@�C������**********************************************************************************************************************
    open (unit=8,form='unformatted',file='phys.q')
    read (8) nblocks
    allocate (Blk_phys(nblocks))
    read (8) (Blk_phys(m)%max(1),Blk_phys(m)%max(2),Blk_phys(m)%max(3),m=1,nblocks)

    do m = 1,nblocks
        allocate (Blk_phys(m)%q(Blk_phys(m)%max(1),Blk_phys(m)%max(2),Blk_phys(m)%max(3),5))
        read (8) mach,alpha,re,t
        read (8) &
            ((((Blk_phys(m)%q(i,j,k,n),i=1,Blk_phys(m)%max(1)),j=1,Blk_phys(m)%max(2)),k=1,Blk_phys(m)%max(3)),n=1,5)
    end do
    close (8)

! �L������**********************************************************************************************************************
    do m = 1,nblocks
        do k = 1,Blk_phys(m)%max(3)
            do j = 1,Blk_phys(m)%max(2)
                do i = 1,Blk_phys(m)%max(1)
                    Blk_phys(m)%q(i,j,k,1) = Blk_phys(m)%q(i,j,k,1)*rho0
                    Blk_phys(m)%q(i,j,k,2) = Blk_phys(m)%q(i,j,k,2)*rho0*c0
                    Blk_phys(m)%q(i,j,k,3) = Blk_phys(m)%q(i,j,k,3)*rho0*c0
                    Blk_phys(m)%q(i,j,k,4) = Blk_phys(m)%q(i,j,k,4)*rho0*c0
                    Blk_phys(m)%q(i,j,k,5) = Blk_phys(m)%q(i,j,k,5)*rho0*c0*c0
                end do
            end do
        end do
    end do

! �t�@�C���o��**********************************************************************************************************************

    open (unit=8,form='unformatted',file='phys_Paraview.q')
    write (8) nblocks
    write (8) (Blk_phys(m)%max(1),Blk_phys(m)%max(2),Blk_phys(m)%max(3),m=1,nblocks)

    do m = 1,nblocks
        write (8) mach,alpha,re,t
        write (8) &
            ((((Blk_phys(m)%q(i,j,k,n),i=1,Blk_phys(m)%max(1)),j=1,Blk_phys(m)%max(2)),k=1,Blk_phys(m)%max(3)),n=1,5)
    end do
    close (8)
    deallocate (Blk_phys)
! �����I�� *************************************************************************************************************************
end subroutine ConvertPhys

! ��`�I�� *************************************************************************************************************************

end program Main
