subroutine decompose_string(string,array)
character*32 array(*)
real tempd
character (*) string
character *256 string2
character *32 temp
integer a1,b1,c1
integer j
j=1

string2=string

do while (len_trim(string2)>0)

if (string2(1:1)==" ") then
string2=string2(2:len_trim(string2))
else
do i=1,len(string2)
    if (string2(i:i)==" ") then
        temp=string2(1:i)
        a1=len_trim(string2)
        string2=string2(i+1:a1)
        array(j)=temp
        array(J+1)=""

        j=j+1
        exit 
    end if
end do
end if
end do
end subroutine