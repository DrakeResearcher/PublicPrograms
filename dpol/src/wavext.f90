!? 28/6/2023 Designed by Eric and Evan

module wavExt

contains

    !? abstraction on top of open() designed for ease of use.
    !! assumes all POL files are in a folder in the working directory called "finite-mass"
    !! assumes all POW files are in a folder in the working directory called "infinite-mass"
    
    !* designed to make migrating to new directory structure as quick as possible:
    !* to convert old code to use this: (assuming your file descriptor is 1)
    !* (replace "open" with "call routedOpen" and youre done!) for example:
    !* open(1, file="FILENAME", status="STATUS)
    !*       becomes
    !* call routedOpen(1, file="FILENAME", status="STATUS)

    subroutine routedOpen(fd, file, form, status)

        integer :: fd, size, foldersize, fullsize, period     
        character :: file*(*), status*(*)
        character :: ext*3
        character(:), allocatable :: path

        character, optional :: form*(*)

        logical :: truncated, tripled


        size = len(file)
        period = findDot(file)
        
        ext = file(period+1:period+3)

        
        if (ext .eq. 'POL') then ! finite mass
            call getSpecs(file, truncated, tripled)

            if (tripled) then
                if (truncated) then
                    foldersize = len("../wave/finite-mass/tripled/truncated/")
                    fullsize = foldersize + size
        
                    allocate(character(fullsize) :: path)
                    path(:foldersize) = "../wave/finite-mass/tripled/truncated/"
                else
                    foldersize = len("../wave/finite-mass/tripled/not-truncated/")
                    fullsize = foldersize + size
        
                    allocate(character(fullsize) :: path)
                    path(:foldersize) = "../wave/finite-mass/tripled/not-truncated/"
                end if
            else
                if (truncated) then
                    foldersize = len("../wave/finite-mass/doubled/truncated/")
                    fullsize = foldersize + size
        
                    allocate(character(fullsize) :: path)
                    path(:foldersize) = "../wave/finite-mass/doubled/truncated/"
                else
                    foldersize = len("../wave/finite-mass/doubled/not-truncated/")
                    fullsize = foldersize + size
        
                    allocate(character(fullsize) :: path)
                    path(:foldersize) = "../wave/finite-mass/doubled/not-truncated/"
                end if
            end if       

        else if (ext .eq. "POW") then ! infitite mass
            call getSpecs(file, truncated, tripled)

            if (tripled) then
                if (truncated) then
                    foldersize = len("../wave/infinite-mass/tripled/truncated/")
                    fullsize = foldersize + size
        
                    allocate(character(fullsize) :: path)
                    path(:foldersize) = "../wave/infinite-mass/tripled/truncated/"
                else
                    foldersize = len("../wave/infinite-mass/tripled/not-truncated/")
                    fullsize = foldersize + size
        
                    allocate(character(fullsize) :: path)
                    path(:foldersize) = "../wave/infinite-mass/tripled/not-truncated/"
                end if
            else
                if (truncated) then
                    foldersize = len("../wave/infinite-mass/doubled/truncated/")
                    fullsize = foldersize + size
        
                    allocate(character(fullsize) :: path)
                    path(:foldersize) = "../wave/infinite-mass/doubled/truncated/"
                else
                    foldersize = len("../wave/infinite-mass/doubled/not-truncated/")
                    fullsize = foldersize + size
        
                    allocate(character(fullsize) :: path)
                    path(:foldersize) = "../wave/infinite-mass/doubled/not-truncated/"
                end if
            end if

        else if (ext .eq. "DAT" .or. ext .eq. "dat") then

            foldersize = len("data/")
            fullsize = foldersize + size

            allocate(character(fullsize) :: path)
            path(:foldersize) = "data/"

        else if (ext .eq. "OUT" .or. ext .eq. "out") then

            foldersize = len("out/")
            fullsize = foldersize + size

            allocate(character(fullsize) :: path)
            path(:foldersize) = "out/"

        else 

            foldersize = len("data/other/")
            fullsize = foldersize + size

            allocate(character(fullsize) :: path)
            path(:foldersize) = "data/other/"

        end if

        path = path(:foldersize) // file
        print *, "==> Looking for file: " // path

        if(present(form)) then    
            open(fd, file=path, form=form, status=status)
        else 
            open(fd, file=path, status=status)
        end if

    end subroutine

    function findDot(str)

        integer findDot
        character str*(*)

        do i=1, len(str)
            if (str(i:i) .eq. ".") then
                findDot = i
                return
            end if
        end do

        findDot = 1

    end function

    subroutine getSpecs(filename, truncated, tripled)

        character :: filename*(*)
        character :: file_in_file*12, line*200, trunc_line*32
        logical :: truncated, tripled

        integer :: first_num, trunc_num, reason

        tripled = .false.
        truncated = .false.

        open(1, file='data/dall2020.dat', status='unknown')
        
        read(1, '(I1)') first_num

        if (first_num .eq. 3) tripled = .true.
        
        print *, filename

        do while (.true.)
            read(1, "(A32)", iostat=reason) line
            if (reason.ne.0) exit

            file_in_file = line(11:22)

            if (file_in_file.eq.filename .or. file_in_file(:11).eq.filename) then
                exit
            end if

        end do

        read(1, "(A24)")
        read(1, "(A24)") trunc_line
        
        trunc_line = trunc_line(19:20)
        read(trunc_line, '(I2)') trunc_num
        
        if (trunc_num.gt.0) truncated = .true.
        close(1)
    end subroutine

end module
