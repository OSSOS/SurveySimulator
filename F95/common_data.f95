module common_data
implicit none
logical, save :: first
integer, save :: iff
character(len=1024), save :: last_surnam
    data first /.true./
    data iff /0/
    data last_surnam /' '/
end module common_data

