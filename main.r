#AUTHOR: VISWAMBHAR REDDY YASA
#MATRICULATION NUMBER :65074
#STOCHASTIC FOR MATERIAL SCIENCE : PROGRAMMING ASSIGNMENT

cat(' PROGRAMMING ASSIGNMENT \n')
cat('STOCHASTIC FOR MATERIAL SCIENCE\n')
cat('NAME : VISWAMBHAR REDDY YASA \n')
cat(' MATRICULATION NUMBER :65074 \n')

cat('select any one of them : task1, task2,task3,task4 \n')
task<-readline(prompt = "Enter the task :")
result = switch(  
  task,  
  "task1"= source('task1_b.r'),  
  "task2"= source('task2.r'),  
  "task3"= source('task3.r'),    
  "task4"= source('task4.r'),  
)
