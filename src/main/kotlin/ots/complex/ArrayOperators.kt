package ots.complex
/**
 * Класс для работы с массивами комплексеных чисел. Комплекснное число из класса Complex
 */

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
 * Сумма действиетельных частей массива
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
    val out: Array<Double> = Array(size) { 0.0 }
    var i=0
    forEach {
        out[i] += it.mod
        i++
    }
    return out
}


/**
 * Функция возвращает обратную матрицу к заданной по методу Гаусса-Жордана
 * @param mat0 входная комплексная матрица
 * @return обратная комплексная матрица
 */
fun findInvMatGaussJordan(mat0: Array<Array<Complex>>): Array<Array<Complex>> {
    val order = mat0.size
    val mat = Array(order) { Array(2 * order) { ZERO } }

    // Заполняем продолженную матрицу: левая часть = mat0, правая часть = единичная матрица
    for (i in 0 until order) {
        for (j in 0 until order) {
            mat[i][j] = mat0[i][j]
        }
        mat[i][i + order] = ONE
    }
    // Прямой ход метода Гаусса-Жордана
    for (i in 0 until order) {
        // Поиск строки с максимальным элементом в текущем столбце
        var maxRow = i
        for (k in i + 1 until order) {
            if (mat[k][i].mod > mat[maxRow][i].mod) {
                maxRow = k
            }
        }
        // Перестановка строк, если необходимо
        if (maxRow != i) {
            val temp = mat[i]
            mat[i] = mat[maxRow]
            mat[maxRow] = temp
        }
        // Проверка на вырожденность матрицы
        if (mat[i][i].mod < 1e-10) {
            return emptyArray()
        }
        // Нормализация текущей строки
        val pivot = mat[i][i]
        for (j in i until 2 * order) {
            mat[i][j] = mat[i][j] / pivot
        }
        // Обнуление других элементов в текущем столбце
        for (k in 0 until order) {
            if (k != i) {
                val factor = mat[k][i]
                for (j in i until 2 * order) {
                    mat[k][j] = mat[k][j] - mat[i][j] * factor
                }
            }
        }
    }
    // Извлечение обратной матрицы из правой части
    val matInv = Array(order) { Array(order) { ZERO } }
    for (i in 0 until order) {
        for (j in 0 until order) {
            matInv[i][j] = mat[i][j + order]
        }
    }
    return matInv
}

/**
 * Функция реализует матричное умножение двух комплексных матриц
 * @param matrixA - комплексная матрица
 * @param matrixB - комплексная матрица
 * @return обратная комплексная матрица
 */
fun multiplyComplexMatrices(
    matrixA: Array<Array<Complex>>,
    matrixB: Array<Array<Complex>>
): Array<Array<Complex>> {
    // Проверка совместимости размеров матриц
    val aRows = matrixA.size
    val aCols = matrixA[0].size
    val bRows = matrixB.size
    val bCols = matrixB[0].size

    if (aCols != bRows) return emptyArray()
    // Создание результирующей матрицы
    val result = Array(aRows) { Array(bCols) { ZERO } }
    // Вычисление произведения матриц
    for (i in 0 until aRows) {
        for (j in 0 until bCols) {
            var sum = ZERO
            for (k in 0 until aCols) {
                sum += matrixA[i][k] * matrixB[k][j]
            }
            result[i][j] = sum
        }
    }
    return result
}