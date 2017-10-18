   program convert
  
!***   nno            : number of nodes
!***   nrtri          : number of elements
!***   nseg           : number of boundary elemenets
!***   itri(1..4,m)   : (1,m) subdomain index for element m (not used here)
!***                    (2..4,m) global number of vertex 1,2,3 for element m
!***   vois(1..3,m)   : element number of triangle opposite to vertex 1..3
!***                    in element m
!***   ifront(1..4,m) : (1..2,m) global node number of vertex (1..2) for
!***                             boundary element m
!***                    (3,m)    number of element to which boundary
!***                             element m belongs
!***                    (4,m)    boundary index for boundary element m
!***                             (not used here)
!***   coeg(1..2,n)   : x,y coordinates of node n
!***   
!***
!***   Convention : the outer boundary is described clockwise, inner
!***                boundaries are described counterclockwise
!***
!***   nno, nrtri and nseg are stored in commonblock.inc
!***
!***   the preset values maxnod, maxele and maxseg in tables.inc 
!***   are checked


 
   implicit none
   integer :: nno, nw, nrtri, nws, nseg, DOF
   integer :: i,j,tecfile, m, n, metisfile, metisfile_dual, progcfile
   integer, dimension(:,:), allocatable  ::itri,ifront,vois
   real*8, dimension(:,:), allocatable :: coeg

   open(1, file = 'mesh.data', status='old', form = 'formatted')
   read (1,'(1x,5(i7,1x))') nno, nw, nrtri, nws, nseg
      
   allocate(itri(4,nrtri))
   allocate(vois(3,nrtri))
   allocate(ifront(4,nseg)) 
   allocate(coeg(2,nno))

   do m = 1, nrtri
      read(1, '(1x,7(i7,1x))')  (itri(i,m), i = 2,4),&
     &                              (vois(i,m), i = 1,3), itri(1,m)
   enddo
   do m = 1, nseg
       read(1, '(1x,4(i7,1x))')   (ifront(i,m), i = 1,4)
   enddo

   do m = 1, nno
       read(1, *)  (coeg(i,m), i = 1,2)
   enddo

   close(1)     




   ! write the mesh in format tecplot: mesh.plt

   tecfile = 10
   open(tecfile, file = 'file.plt', status='unknown', form = 'formatted')
   write(tecfile,*) 'TITLE      =  "2D FE Unstructured grid data" '
   write(tecfile,*) 'VARIABLES = "X", "Y"'
   write(tecfile,*) 'ZONE T="P1", F=FEPOINT, N=',nno,' E=',nrtri,' ET=TRIANGLE'
   
   
   !write the nodes coordinates
   do n= 1,nno
      write(tecfile,'( 2(ES16.8) )') coeg(1,n),coeg(2,n)
   enddo
   
   do n=1, nrtri
      write(tecfile,*) itri(2,n), itri(3,n), itri(4,n)
   enddo
   
   close(tecfile)
      
   ! write the informations needed by metis to create the graph: formetis.data
   metisfile_dual = 10
   open(metisfile_dual,file='dualformetis.dat',status='unknown',form='formatted')
   write(metisfile_dual,*) nrtri, 1 !number of triangles, type of element: 1:triangle,2:tetrahedra,3:hexahedra,4:quadrilateral

   do n=1,nrtri
      write(metisfile_dual,*) itri(2:4,n) !write the element connectivity table (for a mesh partition based on elements : 
                                          !for the creation of the dual graph )
   enddo   

   close(metisfile_dual)

   ! write mesh informations for PROG_c program: Mesh.Data
   progcfile = 11
   open(progcfile,file='meshprogc.data', status= 'unknown', form = 'formatted')

   write(progcfile,*)nno
   do n=1,nno
      DOF = 0!Interior node
      do m = 1,nseg
         if ( (ifront(1,m) == n) .or. (ifront(2,m) == n) )then
              DOF = -1 !Dirichlet Boundary condition
         endif
      enddo
      write(progcfile,*) DOF,coeg(1,n),coeg(2,n)
   enddo
   write(progcfile,*) nrtri
   do m=1,nrtri
      write(progcfile,*) itri(2:4,m) 
   enddo
   close(progcfile)

   deallocate(itri);deallocate(coeg);deallocate(ifront);deallocate(vois);

end program convert
