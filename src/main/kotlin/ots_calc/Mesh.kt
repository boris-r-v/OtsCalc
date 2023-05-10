package ots_calc
import org.kotlinmath.*
//typealias Real = Complex

/**
* Класс сетки
*
* @param X_beg начальная координата участка расчёта, км
* @param X_end конечная координата участка расчёта, км
* @param dX шаг сетки, км
* @property mesh_N количество узлов сетки
* @property X массив узлов сетки с привязкой к координатам
*/
class Mesh(
    val X_beg: Double,
    val X_end: Double,
    val dX: Double,
)
{
    val mesh_N: Int
    private val X: DoubleArray
    init {
        mesh_N = ((X_end - X_beg) / dX).toInt() + 1
        println ("Mesh $mesh_N")
        X = DoubleArray(mesh_N)
        X[0] = X_beg
        for (i in 1 until mesh_N)
            X[i] = X[i - 1] + dX
    }
    /**
     * Возвращает количесво услов сетки
     */
    fun size(): Int{
        return mesh_N
    }
    /**
     * Возвращает элемент сетки по указанному индексу
     */
    fun get(idx: Int): Double
    {
        return X[idx]
    }
    /**
    * Возвращает элементы сетки по индексам укааным во входном массиве
    */
    fun get(idx: IntArray): DoubleArray
    {
        return doubleArrayOf(X[idx[0]], X[idx[1]])
    }
    /**
     * находит ближайший индекс элемента сетки по заданной координате  X
     */
    fun find_near_index_over_mesh(X: Double): Int {
        /* если выходит за левую границу возвращает -1*/
        if (X < this.X[0]) {
            throw trackOutOfMeshException("trackOutOfMeshException point $X out of left edge ${this.X[0]} of mesh")
        }
        /* если выходит за правую границу возвращает -2 */
        if (X > this.X[this.X.size - 1]) {
            throw trackOutOfMeshException("trackOutOfMeshException point $X out of right edge ${this.X[this.X.size - 1]} of mesh")
        }

        return Math.round((X - this.X[0]) / dX).toInt()
    }
    /**
    * находит два ближайших индекса элемента сетки по заданной координате  X
    *          возвращает int[2] - номер ближ элемента слева и справа
    */
    fun find_2near_index_over_mesh(X: Double ): IntArray {
        if (X < this.X[0]) { // если выходит за левую границу
            throw trackOutOfMeshException("trackOutOfMeshException point $X out of left edge ${this.X[0]} of mesh")
        }
        if (X > this.X[this.X.size - 1]) { // если выходит за правую границу
            throw trackOutOfMeshException("trackOutOfMeshException point $X out of right edge ${this.X[this.X.size - 1]} of mesh")
        }
        return intArrayOf(  Math.floor((X - this.X[0]) / dX).toInt(),
                            Math.ceil((X - this.X[0]) / dX).toInt() )
    }
    /**
     * метод распределяет функцию заданную таблицей table[][] на сетку
     *  table[][] это таблица в каждой строке которой содержится два элемента координата км и значение функции
     *  координата в каждой следующей строке должна возрастать
     *  table[][] должен содержать как минимум одну строку,
     *  последняя строка считается до конца сетки не зависимо от координаты в ней
     */
    private fun distribute_function_over_mesh(table: Array<PV>): Array<Real> {
        val N = X.size
        val out: Array<Real> = Array( N ){0.R}
        if (table.size == 1) { // если в таблице только одна строка, то по всей сетке одно значение
            for (i in 0 until N) {
                out[i] = table[0].value
            }
            return out
        }
        //если как минимум две строки
        val N_index = table.size - 1
        var k = 0
        val indexes =
            IntArray(N_index) // создадим массив индексов по сетке для каждой строки за исключением последней
        for (i in 0 until N_index) {
            indexes[i] = find_near_index_over_mesh(table[i].point)
            if (indexes[i] >= 0) {
                while (k <= indexes[i]) {
                    out[k] = table[i].value
                    k++
                }
            }
        }
        while (k < N) {
            out[k] = table[N_index].value
            k++
        }
        return out
    }
    /**
     * создание 3x-Диганальной матрицы в ленточном виде
     * это матрица  состоит из трех строк: нижняя диагональ, главная диагональ и верхняя диагональ
     * на вход принимает экземпляр вложенного класса Track
     */
    fun create_3diag_matrix_band(track: Track): Array<Array<Real>> {
        val N = mesh_N
        var index: Int
        // Распределение по сетке функций сопротивления рельса Ом/км и переходное сопротивление рельс-земля Ом*км
        val r = distribute_function_over_mesh(track.r)
        val rp = distribute_function_over_mesh(track.rp)
        // три вектора: нижняя (под главной) диагональ, главная диагональ, верхняя (над главной диагональю). Эти вектора будут иметь размерность проводимости См
        val diag_dw = Array<Real>(N){0.R}
        val diag = Array<Real>(N){0.R}
        val diag_up = Array<Real>(N){0.R}
        // шаг по сетке [км]
        val dX = dX
        // волновое сопротивление Rv0 [Ом] и сопротивление r0 [Ом] элемента рельса в начале
        val r0: Real
        val rn_mn: Real
        // тоже самое в конце
        var r_i_mn: Real
        var r_i: Real
        var rp_i: Real // сопротивление i-ого элемента рельса и его сопротивление заземления все в Ом

        //формируем массив точечных сосредоточенных сопротивлений в рельсах
        val R_tch = Array<Real>(N - 1){0.R}
        for (i in track.Rtch.indices) {
            index = find_2near_index_over_mesh(track.Rtch[i].point)[0]
            R_tch[index] += track.Rtch[i].value
        }
        // далее все распределенные параметры дискретезуются по сетке в сосредоточенные для каждого элемента
        //из Ом/км -> Ом (вдоль пути); Ом*км -> Ом (на землю)
        if (track.Rv0 == -1.0.R) { // если волновое сопротивление =-1, то определим автоматически
            track.Rv0 = sqrt(r[0] * rp[0]) // sqrt(Ом/км*Ом*км)=sqrt(Ом*Ом)=Ом
        }
        if (track.Rvn == -1.0.R) {
            track.Rvn = sqrt(r[N - 1] * rp[N - 1]) // sqrt(Ом/км*Ом*км)=sqrt(Ом*Ом)=Ом
        }
        r0 = 0.5 * (r[0] + r[1]) * dX + R_tch[0] // Ом/км*км=Ом
        rn_mn = 0.5 * (r[N - 2] + r[N - 1]) * dX + R_tch[N - 2]

        //заполняем первый столбец трех диагоналей
        //все элементы трехдиагональной матрицы имеют размерность См=1/Ом
        diag_up[0] = 0.0.R
        diag[0] = 1 / track.Rv0 + 1 / r0 + 0.5 * dX / rp[0] //  1/Ом+1/Ом+км/(Ом*км)=См
        diag_dw[0] = -1 / r0 //  1/Ом
        //заполняем последний столбец трех диагоналей
        diag_up[N - 1] = -1 / rn_mn //  1/Ом=См
        diag[N - 1] = 1 / track.Rvn + 1 / rn_mn + 0.5 * dX / rp[N - 1] // 1/Ом+1/Ом+км/(Ом*км)=См
        diag_dw[N - 1] = 0.0.R
        //заполянем остальные столбцы
        for (i in 1 until N - 1) {
            r_i_mn = 0.5 * (r[i] + r[i - 1]) * dX + R_tch[i - 1] // Ом/км*км=Ом
            r_i = 0.5 * (r[i] + r[i + 1]) * dX + R_tch[i] // Ом/км*км=Ом
            rp_i = rp[i] / dX // Ом*км/км=Ом
            diag[i] = 1 / r_i_mn + 1 / r_i + 1 / rp_i //заполняем главную диагональ 1/Ом=См
            diag_dw[i - 1] = -1 / r_i_mn //заполняем нижнюю диагональ
            diag_up[i + 1] = -1 / r_i //заполняем верхнюю диагональ
        }
        diag_up[1] = -1 / r0 // дозаполняем нижнюю и верхнюю диагональ строка 2 и предпоследняя
        diag_dw[N - 2] = -1 / rn_mn

        // добавляем проводимости заземлителей zaz в главную диагональ
        for (i in track.zaz.indices) {
            index = find_near_index_over_mesh(track.zaz[i].point)
            diag[index] += 1 / track.zaz[i].value
        }
        return arrayOf(diag_dw, diag, diag_up)
    }
}