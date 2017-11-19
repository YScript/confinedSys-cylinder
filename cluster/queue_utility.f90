MODULE QUEUE_UTILITY
USE global_para
IMPLICIT NONE
PRIVATE
INTEGER,PARAMETER ::  TOP = 1000000
INTEGER, SAVE :: current = 0
INTEGER, SAVE :: que(TOP)
PUBLIC que_clear, que_push, que_pop, que_empty, que_front, que_size

CONTAINS

subroutine que_clear   !��ն���
current = 0
end subroutine que_clear
!*****************************


subroutine que_push(value) !ֻҪ����que_push�ӳ�����ô��ǰ��������+1����ǰ���ݶ�Ӧֵ��value
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



subroutine que_pop !ֻҪ����que_pop���ͻ�ѵ�һ��ֵ�Ӷ�����ȡ������������Ŀcurrent��1�����Ұ��¸�ֵ������һ����que(1)=que(2)��que(2)=que(3).......
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



integer function que_empty() ! �ж϶����Ƿ�Ϊ�գ������ǿյ���ֵΪ1�����в��ǿյ���ֵΪ0 
if (current<=0)then           
que_empty = 1
else
que_empty = 0                 
endif
end function que_empty
!*****************************




integer function que_front()  !�Ѷ��������Ƚ�ȥ��ֵҲ���ǵ�һ��ֵque(1)����que_front��Ҳ�����Ȱ�que(1)���ߣ�����que_front���н������ļ��㡣���У��Ƚ��ȳ�����ջ���Ƚ������
if(current<=0)then            
write(*,*)'que empty'         
stop
return
endif
que_front = que(1)             
end function que_front
!******************************




integer function que_size()    !���еĴ�СҲ���Ƕ��е�����
que_size = current
end function que_size 
!******************************


END MODULE QUEUE_UTILITY
