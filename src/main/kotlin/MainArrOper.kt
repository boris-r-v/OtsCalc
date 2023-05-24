import ots.complex.*
import java.util.*

fun main(args: Array<String>) {
    val ar1 = arrayOf( 2.R, 3.R)
    val ar2 = arrayOf( 20.R, 10.R )

    println("ar1=${Arrays.deepToString(ar1)} ar2=${Arrays.deepToString(ar2)} ar1+ar2=${Arrays.deepToString(ar1+ar2)} ")
    println("ar1=${Arrays.deepToString(ar1)} ar2=${Arrays.deepToString(ar2)} ar1*ar2=${Arrays.deepToString(ar1*ar2)} ")
    println("ar1=${Arrays.deepToString(ar1)} ar2=${Arrays.deepToString(ar2)} ar1-ar2=${Arrays.deepToString(ar1-ar2)} ")
    println("ar1=${Arrays.deepToString(ar1)} ar2=${Arrays.deepToString(ar2)} ar1/ar2=${Arrays.deepToString(ar1/ar2)} ")
    println("ar1=${Arrays.deepToString(ar1)}  ar1*2.0.R=${Arrays.deepToString(ar1*2.0.R)} ")


}