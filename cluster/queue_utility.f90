MODULE QUEUE_UTILITY
USE global_para
IMPLICIT NONE
PRIVATE
INTEGER,PARAMETER ::  TOP = 1000000
INTEGER, SAVE :: current = 0
INTEGER, SAVE :: que(TOP)
PUBLIC que_clear, que_push, que_pop, que_empty, que_front, que_size

CONTAINS

subroutine que_clear   !清空队列
current = 0
end subroutine que_clear
!*****************************


subroutine que_push(value) !只要引用que_push子程序，那么当前数据总数+1，当前数据对应值是value
integer value
if(current>TOP)then
write(*,*)'que full!'
stop
return
endif
current = current + 1
que(current) = value
end subroutine que_push  
!**************************



subroutine que_pop !只要引用que_pop，就会把第一个值从队列中取出，队列总数目current—1，并且把下个值赋给上一个：que(1)=que(2)，que(2)=que(3).......
integer value, i
if(current<=0)then
write(*,*)'que empty!'
stop
return
endif
value = que(1)        
if (current>1)then
  do i = 1, current-1
   que(i) = que(i+1)
  enddo
endif
current = current - 1    
end subroutine que_pop
!*****************************



integer function que_empty() ! 判断队列是否为空，队列是空的则赋值为1，队列不是空的则赋值为0 
if (current<=0)then           
que_empty = 1
else
que_empty = 0                 
endif
end function que_empty
!*****************************




integer function que_front()  !把队列中最先进去的值也就是第一个值que(1)赋给que_front，也就是先把que(1)拿走，交给que_front进行接下来的计算。队列：先进先出。堆栈：先进后出。
if(current<=0)then            
write(*,*)'que empty'         
stop
return
endif
que_front = que(1)             
end function que_front
!******************************




integer function que_size()    !队列的大小也就是队列的容量
que_size = current
end function que_size 
!******************************


END MODULE QUEUE_UTILITY
