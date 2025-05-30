def fun1(): # outer function
    a = 45
    
    def fun2(): # inner function
       # nonlocal a  # allows modification of `a` from fun1
        a=54
        print(a)
    
    fun2()
    print(a)

fun1()
