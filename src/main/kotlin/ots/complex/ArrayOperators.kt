package ots.complex

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