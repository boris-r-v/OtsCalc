package ots.calc
import ots.complex.*
import kotlin.math.ceil
import kotlin.math.floor
import kotlin.math.roundToInt

//typealias Real = Complex

/**
* Класс сетки
*
* @param startX начальная координата участка расчёта, км
* @param endX конечная координата участка расчёта, км
* @param dX шаг сетки, км
* @property meshN количество узлов сетки
* @property x массив узлов сетки с привязкой к координатам
*/
class Mesh(
    val startX: Double,
    val endX: Double,
    val dX: Double,
)
{
    val meshN: Int = ((endX - startX) / dX).toInt() + 1
    internal val x: DoubleArray = DoubleArray(meshN)

    init {
        x[0] = startX
        for (i in 1 until meshN)
            x[i] = x[i - 1] + dX
    }
    /**
     * Возвращает количесво услов сетки
     */
    fun size(): Int{
        return meshN
    }
    /**
     * Возвращает элемент сетки по указанному индексу
     */
    fun get(idx: Int): Double
    {
        return x[idx]
    }
    /**
    * Возвращает элементы сетки по индексам укааным во входном массиве
    */
    fun get(idx: IntArray): DoubleArray
    {
        return doubleArrayOf(x[idx[0]], x[idx[1]])
    }
    /**
     * находит ближайший индекс элемента сетки по заданной координате  X
     */
    fun findNearIndexOverMesh(x: Double): Int {
        /* если выходит за левую границу возвращает -1*/
        if (x < this.x[0]) {
            throw TrackOutOfMeshException("trackOutOfMeshException point $x out of left edge ${this.x[0]} of mesh")
        }
        /* если выходит за правую границу возвращает -2 */
        if (x > this.x[this.x.size - 1]) {
            throw TrackOutOfMeshException("trackOutOfMeshException point $x out of right edge ${this.x[this.x.size - 1]} of mesh")
        }

        return ((x - this.x[0]) / dX).roundToInt()
    }
    /**
    * находит два ближайших индекса элемента сетки по заданной координате  X
    *          возвращает int[2] - номер ближ элемента слева и справа
    */
    fun find2nearIndexOverMesh(x: Double ): IntArray {
        if (x < this.x[0]) { // если выходит за левую границу
            throw TrackOutOfMeshException("trackOutOfMeshException point $x out of left edge ${this.x[0]} of mesh")
        }
        if (x > this.x[this.x.size - 1]) { // если выходит за правую границу
            throw TrackOutOfMeshException("trackOutOfMeshException point $x out of right edge ${this.x[this.x.size - 1]} of mesh")
        }
        return intArrayOf(  floor((x - this.x[0]) / dX).toInt(),
                            ceil((x - this.x[0]) / dX).toInt() )
    }
    /**
     * метод распределяет функцию заданную таблицей table[][] на сетку
     *  table[][] это таблица в каждой строке которой содержится два элемента координата км и значение функции
     *  координата в каждой следующей строке должна возрастать
     *  table[][] должен содержать как минимум одну строку,
     *  последняя строка считается до конца сетки не зависимо от координаты в ней
     */
    private fun distributeFunctionOverMesh(table: Array<PV>): Array<Real> {
        val n = x.size
        val out: Array<Real> = Array( n ){0.R}
        if (table.size == 1) { // если в таблице только одна строка, то по всей сетке одно значение
            for (i in 0 until n) {
                out[i] = table[0].value
            }
            return out
        }
        //если как минимум две строки
        val indexN = table.size - 1
        var k = 0
        val indexes = IntArray(indexN) // создадим массив индексов по сетке для каждой строки за исключением последней
        for (i in 0 until indexN) {
            indexes[i] = findNearIndexOverMesh(table[i].point)
            if (indexes[i] >= 0) {
                while (k <= indexes[i]) {
                    out[k] = table[i].value
                    k++
                }
            }
        }
        while (k < n) {
            out[k] = table[indexN].value
            k++
        }
        return out
    }
    /**
     * создание 3x-Диганальной матрицы в ленточном виде
     * это матрица  состоит из трех строк: нижняя диагональ, главная диагональ и верхняя диагональ
     * на вход принимает экземпляр вложенного класса Track
     */
    fun create3diagMatrixBand(track: Track): Array<Array<Real>> {
        val n = meshN
        var index: Int
        // Распределение по сетке функций сопротивления рельса Ом/км и переходное сопротивление рельс-земля Ом*км
        val r = distributeFunctionOverMesh(track.r)
        val rp = distributeFunctionOverMesh(track.rp)
        // три вектора: нижняя (под главной) диагональ, главная диагональ, верхняя (над главной диагональю). Эти вектора будут иметь размерность проводимости См
        val diagDw = Array(n){0.R}
        val diag = Array(n){0.R}
        val diagUp = Array(n){0.R}
        // шаг по сетке [км]
        val dX = dX
        // волновое сопротивление Rv0 [Ом] и сопротивление r0 [Ом] элемента рельса в начале
        val r0: Real
        val rnMn: Real
        // тоже самое в конце
        var rIMn: Real
        var rI: Real
        var rpI: Real // сопротивление i-ого элемента рельса и его сопротивление заземления все в Ом

        //формируем массив точечных сосредоточенных сопротивлений в рельсах
        val rTch = Array(n - 1){0.R}
        for (i in track.Rtch.indices) {
            index = find2nearIndexOverMesh(track.Rtch[i].point)[0]
            rTch[index] += track.Rtch[i].value
        }
        // далее все распределенные параметры дискретезуются по сетке в сосредоточенные для каждого элемента
        //из Ом/км -> Ом (вдоль пути); Ом*км -> Ом (на землю)
        if (track.Rv0 == (-1.0).R) { // если волновое сопротивление =-1, то определим автоматически
            track.Rv0 = sqrt(r[0] * rp[0]) // sqrt(Ом/км*Ом*км)=sqrt(Ом*Ом)=Ом
        }
        if (track.Rvn == (-1.0).R) {
            track.Rvn = sqrt(r[n - 1] * rp[n - 1]) // sqrt(Ом/км*Ом*км)=sqrt(Ом*Ом)=Ом
        }
        r0 = 0.5 * (r[0] + r[1]) * dX + rTch[0] // Ом/км*км=Ом
        rnMn = 0.5 * (r[n - 2] + r[n - 1]) * dX + rTch[n - 2]

        //заполняем первый столбец трех диагоналей
        //все элементы трехдиагональной матрицы имеют размерность См=1/Ом
        diagUp[0] = 0.0.R
        diag[0] = 1 / track.Rv0 + 1 / r0 + 0.5 * dX / rp[0] //  1/Ом+1/Ом+км/(Ом*км)=См
        diagDw[0] = -1 / r0 //  1/Ом
        //заполняем последний столбец трех диагоналей
        diagUp[n - 1] = -1 / rnMn //  1/Ом=См
        diag[n - 1] = 1 / track.Rvn + 1 / rnMn + 0.5 * dX / rp[n - 1] // 1/Ом+1/Ом+км/(Ом*км)=См
        diagDw[n - 1] = 0.0.R
        //заполянем остальные столбцы
        for (i in 1 until n - 1) {
            rIMn = 0.5 * (r[i] + r[i - 1]) * dX + rTch[i - 1] // Ом/км*км=Ом
            rI = 0.5 * (r[i] + r[i + 1]) * dX + rTch[i] // Ом/км*км=Ом
            rpI = rp[i] / dX // Ом*км/км=Ом
            diag[i] = 1 / rIMn + 1 / rI + 1 / rpI //заполняем главную диагональ 1/Ом=См
            diagDw[i - 1] = -1 / rIMn //заполняем нижнюю диагональ
            diagUp[i + 1] = -1 / rI //заполняем верхнюю диагональ
        }
        diagUp[1] = -1 / r0 // дозаполняем нижнюю и верхнюю диагональ строка 2 и предпоследняя
        diagDw[n - 2] = -1 / rnMn

        // добавляем проводимости заземлителей zaz в главную диагональ
        for (i in track.zaz.indices) {
            index = findNearIndexOverMesh(track.zaz[i].point)
            diag[index] += 1 / track.zaz[i].value
        }
        return arrayOf(diagDw, diag, diagUp)
    }
}