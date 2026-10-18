module common_data
use parameters
implicit none
logical, save :: first
integer, save :: iff
character(len=path_len), save :: last_surnam
    data first /.true./
    data iff /0/
    data last_surnam /' '/
end module common_data

