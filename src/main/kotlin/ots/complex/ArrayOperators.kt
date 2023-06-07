package ots.complex

import kotlin.math.abs

operator fun Array<Complex>.plus(b: Array<Complex>): Array<Complex> {
    if ( size != b.size)
        throw Exception("Arrays size`s different can`t produce array summation")
    return Array(size){ get(it) + b[it] }
}
operator fun Array<Complex>.minus(b: Array<Complex>): Array<Complex> {
    if ( size != b.size)
        throw Exception("Arrays size`s different can`t produce array minus")
    return Array(size){ get(it) - b[it] }
}
operator fun Array<Complex>.times(b: Array<Complex>): Array<Complex> {
    if ( size != b.size)
        throw Exception("Arrays size`s different can`t produce array times")
    return Array(size){ get(it) * b[it] }
}
operator fun Array<Complex>.div(b: Array<Complex>): Array<Complex> {
    if ( size != b.size)
        throw Exception("Arrays size`s different can`t produce array dividing")
    return Array(size){ get(it) / b[it] }
  }

operator fun Array<Complex>.times(b: Complex): Array<Complex> {
    return Array(size){ get(it) * b }
}
operator fun Array<Complex>.plus(b: Complex): Array<Complex> {
    return Array(size){ get(it) + b }
}

/**
 * Сумма модулей элементов массива
 */
fun Array<Complex>.modSum(): Double {
    var out = 0.0
    forEach { out += it.mod }
    return out
}
/**
 * Среднее значение модуля массива, без учета угла
 */
fun Array<Complex>.modAvr(): Double {
    return modSum()/size
}
/**
 * Сумма действиетльных частей массива
 */
fun Array<Complex>.reSum(): Double {
    var out = 0.0
    forEach { out += it.re }
    return out
}
/**
 * Среднее действительных частей
 */
fun Array<Complex>.reAvr(): Double {
    return reSum()/size
}

/**
 * Комплексный массив по модулю как Array<Double>
 */
fun Array<Complex>.mod(): Array<Double>{
    val out: Array<Double> = Array<Double>(size) { 0.0 }
    var i=0
    forEach {
        out[i] += it.mod
        i++
    }
    return out
}