//package org.boristest.OpeOverload
/*
проба работы с классом комплексных чисел
 */
import ots.complex.*

data class Point(val x: Double, val y: Double )

operator fun Point.unaryMinus() = Point(-x, -y)
operator fun Point.unaryPlus() = Point(x, y)
operator fun Point.not() = Point(y, x)
operator fun Point.inc() = Point(x+1, y+1)
operator fun Point.dec() = Point(x-1, y-1)
operator fun Point.plus(b: Point) = Point(x+b.x, y+b.y)
operator fun Point.minus(b: Point) = Point(x-b.x, y-b.y)
operator fun Point.times(b: Point) = Point(x*b.x, y*b.y)
operator fun Point.div(b: Point) = Point(x/b.x, y/b.y)
operator fun Point.rem(b: Point) = Point(x%b.x, y%b.y)
operator fun <T:Int> Point.times(b: T) = Point(x*b, y*b)
operator fun <T:Double> Point.times(b: T) = Point(x*b, y*b)
fun Point.getX() = this.x
fun Point.getY(): Int = TODO("Wait info form master")

  fun main_ref(args: Array<String>) {
    println("Hello World!")


    // Try adding program arguments via Run/Debug configuration.
    // Learn more about running applications: https://www.jetbrains.com/help/idea/running-applications.html.
    println("Program arguments: ${args.joinToString()}")
    val point = Point(1.0,1.0)
    val point2 = Point(2.0,2.0)
    var point3 = Point(2.0,2.0)
    println("plus: ${point*point2}, ${--point3*10} ${(++point3).getX()} ${(point3++).getX()} ")


    val z1 = 4 + 3.I
    val z2 = 4 + 3 * I
    val z3 = complex(3, 4)
    val z4 = "4+3i".toComplex()
    println("z1 $z1")
    println("z2 $z2")
    println("z3 $z3")
    println("z4 $z4")

    println("z1*z2 ${z1*z2}")
    println("z2/z3 ${z2/z3}")
    println("z3+z4 ${z3+z4}")
    println("z4-z1 ${z4-z1}")



    //var numbers: Array<Complex?> = arrayOfNulls(2)//emptyArray()

    var numbers: Array<Complex> = arrayOf(4.R, 4+4.I)//emptyArray()

    for(index in numbers.indices) {
        numbers[index]=numbers[index]*2
    }
    numbers.forEach { println(it/2) }

    val string2DArray: Array<Array<String>> = arrayOf(
        arrayOf("apple", "orange", "avocado", "mango", "banana"),
        arrayOf("_", "!", ":", "?"),
        arrayOf("1", "2", "3", "4", "5", "10"))

    // Print the 2D array
    string2DArray.forEach {
        it.forEach { it -> print("$it, ") }
        println()
    }

    // Access an individual String using index notation.
    println("My favorite fruit is: ${string2DArray[0][2]}")

    val complex2d: Array<Array<Complex>> = emptyArray()

    /*
        println("Hello World!")

        val complex2d: Array<Array<Complex>> = arrayOf(
            arrayOf(4.R, 4+1.I),
            arrayOf(5.I, 4*exp(3.I))
        )
        val dd = SomeComplexClass(complex2d, 5.5)
    */
}