module f2py_api
    public :: my_matmult_api

contains

    subroutine my_matmult_api(a, b, n, c)
        ! This is the public interface to the GPU matmult that
        ! we call from python
        use cudafor ! this has to be hear in order to make the module call
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(in) :: b(:, :)
        integer, intent(in) :: n
        real(8), intent(out) :: c(n, n)
        integer, i,j,k
        integer, nThreads 
        ! --- GPU DEVICE PARAMETERS ---
        real(8), allocatable,device :: d_a(:, :), d_b(:, :), d_c(:, :)
        type(dim3) :: grid, tBlock

        ! Variables for timing
        integer :: start, finish, count_rate
        real(8) :: elapsed_time


        ! --- GPU device grid and block size ---
        nThreads = 32
        grid = dim3(ceiling(real(n) / nThreads), ceiling(real(n) / nThreads),1) ! we ask for more from the GPU using the 'ceiling()'
        tBlock = dim3(nThreads, nThreads,1) ! 32 by 32 is the maximum you can ask for
        print*,grid%x,grid%y
        print *,tBlock%x, tBlock%y
        c = 0.0
        ! do i =1,n
        !     do j =1,N
        !         do k = 1,n
        !             c(i,j) = c(i,j)+a(i,k)*b(k,j)
        !         end do 
        !     end do 
        ! end do
        ! --- Allocation ---
        ! allocate
        call system_clock(start, count_rate) ! get start time

        allocate (d_a(n, n), d_b(n, n), d_c(n, n))
        d_a = a
        d_b = b
        d_c = 0
        !call cuda
        call my_matmult_gpu<<<grid, tBlock>>>(d_a, d_b, d_c)
        ! copy memory back
        c = d_c
        ! print *, "c(1,1) = ", c(1, 1)
        ! convert time to seconds and print
        call system_clock(finish) ! get finish time

        elapsed_time = real(finish - start, 8) / real(count_rate, 8)
        write (*, '(a,f9.4,a)') "gpu mat-mat product", elapsed_time, " seconds!"

    end subroutine my_matmult_api


    attributes(global) subroutine my_matmult_gpu(a, b, c)
        implicit none
        real(8) :: a(:, :)
        real(8) :: b(:, :)
        real(8) :: c(:, :)

        ! --- Working vars ---
        real(8) :: tmp
        integer :: n
        integer :: ii, jj ! thread indices
        integer :: kk !loop index
        integer :: blockSize
        integer :: tx, ty, bx, by
        ! real(8), shared :: aTile(32, 32), bTile(32, 32)
        integer :: ll, mm ! loop indices for tiling
        integer :: rowGlob, colGlob
        n = size(a, dim=1)

        ! At a bare minimum, any thread that can do work has at minimum two identifying numbers:
        ! block index
        ! thread index
        ! There is also blockDim which is how many threads are in a block

        ! Convenience
        tx = threadIdx%x
        ty = threadIdx%y
        bx = blockIdx%x
        by = blockIdx%y
        ! ************************************************
        ! !     Naive implementation of matrix multiplication (coalleced too because its fortran :)
        ! ! ************************************************
        ! ! Ok so now...I've made it so there are 2 id numbers within the thread index so there's %x and %y
        ii = tx + (bx - 1) * blockDim%x
        jj = ty + (by - 1) * blockDim%y
        ! ! The blockIdx%x - 1 is because FORTRAN is 1-based
        if (ii <= n .and. jj <= n) then
            tmp = 0.0
            do kk = 1, n
                ! tmp=0.0
                tmp = tmp + a(ii, kk) * b(kk, jj)
            end do
            c(ii, jj) = tmp
        end if


    end subroutine my_matmult_gpu
     

end module f2py_api
