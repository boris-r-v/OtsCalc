package ots.calc

import java.util.*

/** Класс для расчёта мгновенной схемы ОТС
 * на постоянном токе НЕ однородный рельс
 * -----------------!!!НЕ ИСПОЛЬЗУЕТСЯ на переменном токе!!!-------------------
 */
class OTS_DC_HeterogRail_call_kt internal constructor(// количество главных путей и ответвлений
    private val Ntrc: Int, // количество сеток
    private val Nmsh: Int
) {
    /* В комментариях в данном классе используются следующие обозначения:
  Р - рельсы данного пути, рассматривается как длинная линия
  ЭПС - ЭПС по данному пути
  ФОТ - фидера обратного тока (отсасывающие фидера) подключенные к рельсам данного пути
  ЗАЗ - заземлители к рельсам данного пути (имеет смысл учитывать заземлители с малым сопротивлением заземления менее 10 Ом)
  МПС -  междупутные соединители, условное название – фактически это сопротивление, соединяющее две любые точки каких любо путей (длинных линий)
  МКР – метод конечных разностей, сеточный метод дискретизации дифференциальных уравнений (в данном случае уравнение длинной электрической линии) по пространству (здесь вдоль пути)  

  
  Единицы измерения в классе:
  координаты и длины - км;
  напряжение  - В;
  ток - А;
  сопротивление - Ом;
  проводимость - См;
  погонное сопротивление - Ом/км;
  переходное сопротивление Ом*км;
  */
    //-----------------------------------------------------Головной класс переменные и инициализация-------------------------------------------------------------------------
    // переменные дочерних классов
    var computing_settings: Computing_settings
    var err_and_mes: Errors_and_messages
    var meshes // массив сеток
            : Array<Mesh?>
    var tracks // массивы путей
            : Array<Track?>

    //переменные стандартных типов данных
    var mps =arrayOf(doubleArrayOf()) // двумерный массив в каждой строке 5 элементов:
                                    // 0 - номер пути начальной точки подключения,
                                    // 1 - номер пути конечной точки подключения,
                                    // 2 - координата начальной точки подключения км,
                                    // 3 - координата конечной точки подключения км,
                                    // 4 - сопротивление Ом.
    private lateinit var I_mps: DoubleArray // массив токов междупутных соединителей МПС и токов отходящих ветвей в местах соединения с главными путями
    private lateinit var num_mesh_mps: Array<IntArray>

    // двумерный массивы номеров сеток для МПС (в каждой строке начальная и конечная точка)
    private lateinit var num_track_mps: Array<IntArray>

    // двумерный номеров путей для МПС (в каждой строке начальная и конечная точка)
    private lateinit var index_mps: Array<IntArray> //двумерный массивы  индексов узлов по сетке МПС (в каждой строке начальная и конечная точка)
    private lateinit var a_x_find: Array<DoubleArray> // матрица коэффициентов влияния тока во всех МПС на напряжения во всех МПС Ом. По главной диагонали сами на себя.
    private lateinit var U_const: DoubleArray //  массивы разности напряжений от заданных токов в МПС от начальнйо до конечной точки подключения

    // инициализатор основного класса по умолчанию
    init { // количество путей и сеток
        computing_settings = Computing_settings()
        err_and_mes = Errors_and_messages()
        meshes = arrayOfNulls(Nmsh)
        tracks = arrayOfNulls(Ntrc)
        for (i in 0 until Ntrc) { // заполняем массив путей
            tracks[i] = Track()
        }
        for (i in 0 until Nmsh) { // заполняем массив сеток
            meshes[i] = Mesh()
        }
    }

    //-----------------------------------------------------Подклассы, их методы, свойства и инициализаторы-------------------------------------------------------------------------
    /*совмещенный класс для пути
     */
    inner class Track internal constructor() {
        lateinit var r: Array<DoubleArray>

        //функция погонного сопротивления рельсов вдоль пути. В каждой строке два значения: 1- координата км, 2 -величина погонного сопротивление рельсов Ом/км до данной координаты.
        lateinit var rp: Array<DoubleArray>

        // функция переходного сопротивления рельсы-земля вдоль пути. В каждой строке два значения: 1- координата км, 2 -величина переходного сопротивления Ом*км до данной координаты.
        var fot: Array<DoubleArray>

        // таблица подключенных к данному пути ФОТ. В каждой строке два значения: 1- координата точки подключения км, 2 -величина тока ФОТ А. Положительный ток принят при работе ТП в режиме выпрямителя (когда ток из рельсов затекает в ФОТ)
        var eps: Array<DoubleArray>

        // таблица ЭПС находящихся на данном пути. В каждой строке два значения: 1- координата расположения км, 2 -величина тока ЭПС А. Положительный ток принят при работе ЭПС в режиме тяги или хх (когда ток стекает из колесных пар в рельсы)
        var zaz: Array<DoubleArray>

        // таблица ЗАЗ находящихся на данном пути. В каждой строке два значения: 1- координата точки подключения км, 2 -величина сопротивления заземления Ом.
        var R_tch: Array<DoubleArray>

        // таблица сосредоточенных точечных сопротивлений в рельсах на данном пути. В каждой строке два значения: 1- координата расположения км, 2 -величина сосредоточенного сопротивления Ом. (к примеру: разрыв рельсовой сети, дефектный неизолированный стык с высоким сопротивлением, дефект в средней точке ДТ с повышенным сопротивлением).
        var Rv0: Double

        // волновое сопротивление слева Ом
        var Rvn // волновое сопротивление справа Ом
                : Double
        var num_mesh = 0 // номер сетки для данного пути
        lateinit var m3db: Array<DoubleArray>

        // матрица 3диагональная ленточная См
        lateinit var vector_b: DoubleArray

        // вектор правой части А
        lateinit var u: DoubleArray

        // искомый вектор напряжения в узлах сетки  В
        lateinit var i: DoubleArray

        // вектор тока рельсов в узлах сетки А
        lateinit var i_grnd: DoubleArray // вектор тока земле в узлах сетки от данного пути А

        init {
            fot = arrayOf()
            eps = arrayOf()
            zaz = arrayOf()
            R_tch = arrayOf()
            Rv0 =
                -1.0 // присвоим волновое сопротивление в начале и в конце значение -1, что автоматически означает расчёт по величине r и rp
            Rvn = -1.0
        }
    }

    /** Класс расчётная сетка
     * дискретизация вдоль пути рельсов как электрически длинной линии
     * линейная однородная сетка как основа метода конечных разностей
     */
    inner class Mesh // по умолчанию параметры сетки такие
    internal constructor() {
        var X_beg = 0.0
        var X_end //начальная и конечная координата участка расчёта км
                = 20.0
        var dX //шаг сетки по длине в км
                = 0.1
        lateinit var X: DoubleArray // массив координат узлов сетки
        var mesh_N = 0 // количество узлов сетки

        /*метод создает сетку
        за основу берет данные из полей класса: X_beg,X_end, dX
        */
        fun create_mesh() {
            val X0 = dX * Math.round(X_beg / dX).toInt()
            val Xn = dX * Math.round(X_end / dX).toInt()
            mesh_N = ((Xn - X0) / dX).toInt() + 1
            X = DoubleArray(mesh_N)
            X[0] = X0
            for (i in 1 until mesh_N) {
                X[i] = X[i - 1] + dX
            }
        }

        /*метод распределяет функцию заданную таблицей table[][] на сетку
        table[][] это таблица в каждой строке которой содержится два элемента координата км и значение функции
        координата в каждой следующей строке должна возрастать
         table[][] должен содержать как минимум одну строку, последняя строка считается до конца сетки не зависимо от координаты в ней
        */
        private fun distribute_function_over_mesh(table: Array<DoubleArray>): DoubleArray {
            val N = X.size
            val out = DoubleArray(N)
            if (table.size == 1) { // если в таблице только одна строка, то по всей сетке одно значение
                for (i in 0 until N) {
                    out[i] = table[0][1]
                }
                return out
            }
            //если как минимум две строки
            val N_index = table.size - 1
            var k = 0
            val indexes =
                IntArray(N_index) // создадим массив индексов по сетке для каждой строки за исключением последней
            for (i in 0 until N_index) {
                indexes[i] = find_near_index_over_mesh(table[i][0])
                if (indexes[i] >= 0) {
                    while (k <= indexes[i]) {
                        out[k] = table[i][1]
                        k++
                    }
                }
            }
            while (k < N) {
                out[k] = table[N_index][1]
                k++
            }
            return out
        }

        // находит ближайший индекс элемента сетки по заданной координате  X
        fun find_near_index_over_mesh(X: Double): Int {
            if (X < this.X[0]) { // если выходит за левую границу возвращает -1
                return -1
            }
            return if (X > this.X[this.X.size - 1]) { // если выходит за правую границу возвращает -2
                -2
            } else Math.round((X - this.X[0]) / dX).toInt()
        }

        /* находит два ближайших индекса элемента сетки по заданной координате  X
         возвращает int[2] - номер ближайшего элемента слева и справа
         */
        fun find_2near_index_over_mesh(X: Double): IntArray {
            if (X < this.X[0]) { // если выходит за левую границу
                return intArrayOf(-1, -1)
            }
            return if (X > this.X[this.X.size - 1]) { // если выходит за правую границу
                intArrayOf(-2, -2)
            } else intArrayOf(
                Math.floor((X - this.X[0]) / dX).toInt(),
                Math.ceil((X - this.X[0]) / dX).toInt()
            )
        }

        /*создание 3Диганальной матрицы в ленточном виде
        это матрица  состоит из трех строк: нижняя диагональ, главная диагональ и верхняя диагональ
        на вход принимает экземпляр вложенного класса Track
      */
        fun create_3diag_matrix_band(track: Track?) {
            val N = mesh_N
            var index: Int
            val r = distribute_function_over_mesh(track!!.r)
            val rp = distribute_function_over_mesh(track.rp)
            // Распределение по сетке функций сопротивления рельса Ом/км и переходное сопротивление рельс-земля Ом*км
            val diag_dw = DoubleArray(N)
            val diag = DoubleArray(N)
            val diag_up = DoubleArray(N)
            // три вектора: нижняя (под главной) диагональ, главная диагональ, верхняя (над главной диагональю). Эти вектора будут иметь размерность проводимости См
            val dX = dX
            // шаг по сетке [км]
            val r0: Double
            // волновое сопротивление Rv0 [Ом] и сопротивление r0 [Ом] элемента рельса в начале
            val rn_mn: Double
            // тоже самое в конце
            var r_i_mn: Double
            var r_i: Double
            var rp_i: Double // сопротивление i-ого элемента рельса и его сопротивление заземления все в Ом

            //формируем массив точечных сосредоточенных сопротивлений в рельсах
            val R_tch = DoubleArray(N - 1)
            for (i in track.R_tch.indices) {
                index = find_2near_index_over_mesh(track.R_tch[i][0])[0]
                R_tch[index] += track.R_tch[i][1]
            }
            // далее все распределенные параметры дискретезуются по сетке в сосредоточенные для каждого элемента
            //из Ом/км -> Ом (вдоль пути); Ом*км -> Ом (на землю)
            if (track.Rv0 == -1.0) { // если волновое сопротивление =-1, то определим автоматически
                track.Rv0 = Math.sqrt(r[0] * rp[0]) // sqrt(Ом/км*Ом*км)=sqrt(Ом*Ом)=Ом
            }
            if (track.Rvn == -1.0) {
                track.Rvn = Math.sqrt(r[N - 1] * rp[N - 1]) // sqrt(Ом/км*Ом*км)=sqrt(Ом*Ом)=Ом
            }
            r0 = 0.5 * (r[0] + r[1]) * dX + R_tch[0] // Ом/км*км=Ом
            rn_mn = 0.5 * (r[N - 2] + r[N - 1]) * dX + R_tch[N - 2]

            //заполняем первый столбец трех диагоналей
            //все элементы трехдиагональной матрицы имеют размерность См=1/Ом
            diag_up[0] = 0.0
            diag[0] = 1 / track.Rv0 + 1 / r0 + 0.5 * dX / rp[0] //  1/Ом+1/Ом+км/(Ом*км)=См
            diag_dw[0] = -1 / r0 //  1/Ом
            //заполняем последний столбец трех диагоналей
            diag_up[N - 1] = -1 / rn_mn //  1/Ом=См
            diag[N - 1] = 1 / track.Rvn + 1 / rn_mn + 0.5 * dX / rp[N - 1] // 1/Ом+1/Ом+км/(Ом*км)=См
            diag_dw[N - 1] = 0.0
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
                index = find_near_index_over_mesh(track.zaz[i][0])
                diag[index] += 1 / track.zaz[i][1]
            }
            track.m3db = arrayOf(diag_dw, diag, diag_up)
        }

        /*геттер возвращает массив координат узлов сетки
         */
        fun get_X(): DoubleArray {
            return X
        }
    }

    /** класс вычислительных настроек для электрического расчёта ОТС
     */
    inner class Computing_settings // инициализатор по умолчанию
    {
        var convergence_U // допустимая невязка величины напряжения усреднённая по всем точкам граничного условия, в которых происходит итерационный поиск величины тока втекающего в эти точки.
                = 0.01
        var max_number_iterat // максимальное число итераций при расчёте величины тока, втекающего в граничные точки, по методу Ньютона в цикле.
                = 1000
        var initial_damping_factor // начальное коэффициент демпфирования в методе Ньютона при расчёте величины тока, втекающего в граничные точки.
                = 0.7
        lateinit var current_state_solver: DoubleArray // текущее состояние решателя содержит массив из трёх чисел полученных в конце последнего расчёта: 1 - количество итераций, 2- средняя невязка по напряжению, 3 - коэффициент демпфирования

        //возвращает параметры решаттеля в виде массива
        fun get_current_state_solver(): DoubleArray {
            return current_state_solver
        }
    }

    /** класс для ошибок при вычислениях, вводе данных, и сообщений о них
     */
    inner class verify_data // инициализатор по умолчанию - ошибок нет, расчёт не выполнен
    {
        var data_error = false
        var solver_error = false
        var calc_completed //data_error - ошибка в данных (к примеру, неправильная длина массива), solver_error - шибка решателя (к примеру, не достигнута сходимость, исчерпано число итераций), calc_completed  - признак что расчёт выполнен хотя бы раз
                = false
        var messeg_data_error = ""
        var messeg_solver_error // текстовое сообщение об этих ошибках
                = "Расчёт не выполнен"

        // геттеры для всех свойств
        fun get_data_error(): Boolean {
            return data_error
        }

        fun get_solver_error(): Boolean {
            return solver_error
        }

        fun get_messeg_data_error(): String {
            return messeg_data_error
        }

        fun get_messeg_solver_error(): String {
            return messeg_solver_error
        }

        //обнуление для ошибки данных
        fun reset_data_error() {
            data_error = false
            messeg_data_error = ""
        }

        //обнуление для ошибки решателя
        fun reset_solver_error() {
            solver_error = false
            messeg_solver_error = ""
        }
    }
    /** класс для ошибок при вычислениях, вводе данных, и сообщений о них
     */
    inner class Errors_and_messages // инициализатор по умолчанию - ошибок нет, расчёт не выполнен
    {
        var data_error = false
        var solver_error = false
        var calc_completed = false //data_error - ошибка в данных (к примеру, неправильная длина массива), solver_error - шибка решателя (к примеру, не достигнута сходимость, исчерпано число итераций), calc_completed  - признак что расчёт выполнен хотя бы раз
        var messeg_data_error = ""
        var messeg_solver_error = "Расчёт не выполнен" // текстовое сообщение об этих ошибках
        // геттеры для всех свойств
        fun get_data_error(): Boolean {
            return data_error
        }
        fun get_solver_error(): Boolean {
            return solver_error
        }
        fun get_messeg_data_error(): String {
            return messeg_data_error
        }
        fun get_messeg_solver_error(): String {
            return messeg_solver_error
        }
        fun reset_data_error() {
            data_error = false
            messeg_data_error = ""
        }
        fun reset_solver_error() {
            solver_error = false
            messeg_solver_error = ""
        }
    }
    //-----------------------------------------------------Внутренние процедуры базового класса приватные и публичные (за исключением сеттеров и геттеров)-------------------------------------------------------------------------
    //проверка исходных данных
    private fun verify_data(
        verify_I: Boolean,
        verify_Xfot_mps: Boolean,
        verify_Xeps: Boolean
    ) { // verify_I - условие проверки токов; verify_Xfot_mps - проверка коорд ФОТ, МПС, ЗАЗ; verify_Xeps - проверка коорд ЭПС
        err_and_mes.reset_data_error() // очистка data_error

        // проверка по заданию ФОТ и ЭПС и токам
        var Neps_fot_all = 0
        var Isum_fot = 0.0
        var Isum_eps = 0.0
        val deltaI: Double
        if (verify_I) {
            // проверяем чтобы сумма токов ЭПС равнялась сумме токов ФОТ не превышая погрешности
            // и что вообще ЭПС  и ФОТ заданы
            for (i in 0 until Ntrc) {
                Neps_fot_all += tracks[i]!!.fot.size + tracks[i]!!.eps.size
                for (j in tracks[i]!!.fot.indices) {
                    Isum_fot += tracks[i]!!.fot[j][1]
                }
                for (j in tracks[i]!!.eps.indices) {
                    Isum_eps += tracks[i]!!.eps[j][1]
                }
            }
            if (Neps_fot_all == 0) { // если массивы ФОТ и/или ЭПС не заданы
                err_and_mes.data_error = true
                err_and_mes.messeg_data_error += "Массивы ФОТ и/или ЭПС не заданы. "
            }
            deltaI = Math.abs(Isum_fot - Isum_eps)
            if (deltaI / (Math.max(Math.abs(Isum_fot), Math.abs(Isum_eps)) + 1e-6) > 0.05) {
                err_and_mes.messeg_data_error += " Предупреждение: сумма токов ФОТ и ЭПС расходится более чем на 5 %"
            }
        }

        // проверки по сетке и координатам
        var all_points: Array<Array<DoubleArray>>
        var num_mesh1: Int
        var num_mesh2: Int
        var num_track1: Int
        var num_track2: Int
        if (verify_Xfot_mps) {
            // проверяем чтобы  границы сетки заданы корректно
            for (i in 0 until Nmsh) {
                if (meshes[i]!!.X_end - meshes[i]!!.X_beg <= 2 * meshes[i]!!.dX) { // если ошибка в границах сетки и менее двух узлов
                    err_and_mes.data_error = true
                    err_and_mes.messeg_data_error += "Сетка номер $i: границы сетки заданы не корректно, либо получается менее трёх узлов"
                }
            }
            //проверяем чтобы координаты точек ФОТ, ЗАЗ, R_tch укладывались в границы сетки данного пути
            for (i in 0 until Ntrc) {
                all_points = arrayOf(tracks[i]!!.fot, tracks[i]!!.zaz, tracks[i]!!.R_tch)
                num_mesh1 = tracks[i]!!.num_mesh
                for (j in all_points.indices) {
                    for (k in all_points[j].indices) {
                        if (all_points[j][k][0] > meshes[num_mesh1]!!.X_end || all_points[j][k][0] < meshes[num_mesh1]!!.X_beg) {
                            err_and_mes.data_error = true
                            err_and_mes.messeg_data_error += "Путь номер " + (i + 1) + ": координаты точек ФОТ, ЗАЗ или сосред сопротивлен выходят за границы сетки"
                        }
                    }
                }
            }
            // проверка чтобы точки подключения МПС к путям в пределах сетки
            for (i in mps.indices) {
                num_track1 = mps[i][0].toInt()
                num_track2 = mps[i][1].toInt()
                num_mesh1 = tracks[num_track1]!!.num_mesh
                num_mesh2 = tracks[num_track2]!!.num_mesh
                if (mps[i][2] > meshes[num_mesh1]!!.X_end || mps[i][2] < meshes[num_mesh1]!!.X_beg) {
                    err_and_mes.data_error = true
                    err_and_mes.messeg_data_error += "МПС " + i + ": координата начальной точки подключения к пути " + num_track1 + "выходbт за границы сетки"
                }
                if (mps[i][3] > meshes[num_mesh2]!!.X_end || mps[i][3] < meshes[num_mesh2]!!.X_beg) {
                    err_and_mes.data_error = true
                    err_and_mes.messeg_data_error += "МПС " + i + ": координата конечной точки подключения к пути " + num_track2 + "выходbт за границы сетки"
                }
            }
        }
        if (verify_Xfot_mps) {
            //проверяем чтобы координаты точек ЭПС укладывались в границы сетки данного пути
            for (i in 0 until Ntrc) {
                num_mesh1 = tracks[i]!!.num_mesh
                for (j in tracks[i]!!.eps.indices) {
                    if (tracks[i]!!.eps[j][0] > meshes[num_mesh1]!!.X_end || tracks[i]!!.eps[j][0] < meshes[num_mesh1]!!.X_beg) {
                        err_and_mes.data_error = true
                        err_and_mes.messeg_data_error += "Путь номер " + (i + 1) + ": координаты точек ФОТ, ЗАЗ или сосред сопротивлен выходят за границы сетки"
                    }
                }
            }
        }
    }

    /** Метод для решения СЛАУ с 3диагонал ленточной матрицей
     * методом двойной проходки
     * на входе матрица и вектор правой части
     * на выходе вектор ответов
     */
    fun solve_3diag_band(
        matrix_band: Array<DoubleArray>,
        vector_b: DoubleArray
    ): DoubleArray { // matrix_band – трёхдиагональная ленточная матрица, vector_b - вектор правой части
        val N = vector_b.size
        val v = DoubleArray(N)
        val u = DoubleArray(N)
        val out = DoubleArray(N)
        //прямая проходка
        v[0] = matrix_band[2][1] / -matrix_band[1][0]
        u[0] = -vector_b[0] / -matrix_band[1][0]
        for (i in 1 until N - 1) {
            v[i] = matrix_band[2][i + 1] / (-matrix_band[1][i] - matrix_band[0][i - 1] * v[i - 1])
            u[i] =
                (matrix_band[0][i - 1] * u[i - 1] - vector_b[i]) / (-matrix_band[1][i] - matrix_band[0][i - 1] * v[i - 1])
        }
        v[N - 1] = 0.0
        u[N - 1] =
            (matrix_band[0][N - 2] * u[N - 2] - vector_b[N - 1]) / (-matrix_band[1][N - 1] - matrix_band[0][N - 2] * v[N - 2])
        //обратная проходка
        out[N - 1] = u[N - 1]
        for (i in N - 1 downTo 1) {
            out[i - 1] = out[i] * v[i - 1] + u[i - 1]
        }
        return out
    }

    /** Заполняет массивы номеров путей, сеток и индексов узлов для МПС
     * в каждой строке два значения int: для начальной и конечной точки подключения
     */
    private fun find_index_mps() {
        val N = mps.size
        num_track_mps = Array(N) { IntArray(2) } // инициализируем массивы
        num_mesh_mps = Array(N) { IntArray(2) }
        index_mps = Array(N) { IntArray(2) }
        for (i in 0 until N) { //заполняем массивы
            num_track_mps[i][0] = mps[i][0].toInt()
            num_track_mps[i][1] = mps[i][1].toInt()
            num_mesh_mps[i][0] = tracks[num_track_mps[i][0]]!!.num_mesh
            num_mesh_mps[i][1] = tracks[num_track_mps[i][1]]!!.num_mesh
            index_mps[i][0] = meshes[num_mesh_mps[i][0]]!!.find_near_index_over_mesh(mps[i][2])
            index_mps[i][1] = meshes[num_mesh_mps[i][1]]!!.find_near_index_over_mesh(mps[i][3])
        }
    }

    /** Метод инициализирует ОТС
     */
    fun init_OTS() {
        var num_mesh: Int // номер сетки дял текущего пути
        verify_data(true, true, true) // проверка исходных данных полная
        for (i in 0 until Nmsh) {
            meshes[i]!!.create_mesh() //создаём сетку
        }
        find_index_mps() // инициализируем массивы номеров начальной и конечной точки подключения для  МПС
        U_const = DoubleArray(mps.size) // инициализируем массивы постоянной составляющей напряжения
        for (i in 0 until Ntrc) { // для обоих путей
            num_mesh = tracks[i]!!.num_mesh
            tracks[i]!!.vector_b = DoubleArray(meshes[num_mesh]!!.mesh_N) // инициализируем вектор правой части
            meshes[num_mesh]!!.create_3diag_matrix_band(tracks[i]) // рассчитываем 3диагональную ленточную матрицу
            tracks[i]!!.u = DoubleArray(meshes[num_mesh]!!.mesh_N) // инициализация векторов напряжений в рельсах
            tracks[i]!!.i = DoubleArray(meshes[num_mesh]!!.mesh_N) // токов в рельсах
            tracks[i]!!.i_grnd = DoubleArray(meshes[num_mesh]!!.mesh_N) // токов из рельсов в землю
        }
        eval_a_x_find() // рассчитываем матрицу коэффициентов влияния в МПС
    }

    /* Расчет постоянной составляющей напряжения в точках МПС
     */
    private fun eval_U_const() {
        // Рассчитаем напряжение в узлах каждого пути от постоянных источников
        for (i in 0 until Ntrc) {
            tracks[i]!!.vector_b = addValues_vector_b_1node(
                tracks[i]!!.num_mesh,
                tracks[i]!!.vector_b,
                tracks[i]!!.fot,
                -1.0
            ) // в точках ФОТ ток в одном узле с минусом
            tracks[i]!!.vector_b = addValues_vector_b_2node(
                tracks[i]!!.num_mesh,
                tracks[i]!!.vector_b,
                tracks[i]!!.eps,
                1.0
            ) // в точках ЭПС ток  в двух ближайших узлах
            tracks[i]!!.u = solve_3diag_band(tracks[i]!!.m3db, tracks[i]!!.vector_b) //потенциал в рельсах в узлах сетки
            println("tracks[$i]!!.u, ${Arrays.toString( tracks[i]!!.u) } ")
        }

        // по параметрам МПС (номерам путей и координатам точек начала и конца) запишем в U1_const и U2_const  в узлах с МПС  напряжения в рельсах от постоянных источников
        for (j in mps.indices) {
            U_const[j] =
                tracks[num_track_mps[j][0]]!!.u[index_mps[j][0]] - tracks[num_track_mps[j][1]]!!.u[index_mps[j][1]]
        }
        println("U_const:, ${Arrays.toString(U_const)} ")
    }

    /* Расчет коэффициентов влияния в МПС от их самих
     */
    private fun eval_a_x_find() {
        val N_find = mps.size // количестов поисковых точек
        val Imps = 1000.0
        // уловный ток МПС для определения коэффициентов влияния
        var u1: DoubleArray
        var u2: DoubleArray
        // массивы напряжений в пути начальной и конечной точки подключения создаваемые током МПС
        val a1_x_find = Array(N_find) { DoubleArray(N_find) }
        // временные массивы влияния в МПС в начальной и
        val a2_x_find = Array(N_find) { DoubleArray(N_find) } // в конечной
        a_x_find = Array(N_find) { DoubleArray(N_find) } // инициализируем общий массив влияния
        // Перебераем поисковые точки МПС
        for (i in 0 until N_find) {
            tracks[num_track_mps[i][0]]!!.vector_b[index_mps[i][0]] =
                -Imps // задаём ток в начальной точке подключения в данном пути
            tracks[num_track_mps[i][1]]!!.vector_b[index_mps[i][1]] =
                Imps // задаём ток в конечной точке подключения в данном пути
            u1 = solve_3diag_band(
                tracks[num_track_mps[i][0]]!!.m3db,
                tracks[num_track_mps[i][0]]!!.vector_b
            ) //снимаем напряжение по сетке в пути для начальной точки подключения
            u2 = solve_3diag_band(
                tracks[num_track_mps[i][1]]!!.m3db,
                tracks[num_track_mps[i][1]]!!.vector_b
            ) //снимаем напряжение по сетке в пути для конечной точки подключения
            tracks[num_track_mps[i][0]]!!.vector_b[index_mps[i][0]] =
                0.0 // обнуляем ток в текущей начальной точке подключения
            tracks[num_track_mps[i][1]]!!.vector_b[index_mps[i][1]] =
                0.0 // обнуляем ток в текущей конечной точке подключения
            for (j in 0 until N_find) { //снова проходим по всем точкам МПС
                if (num_track_mps[i][0] == num_track_mps[j][0]) { // если номер пути начальной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1_x_find[i][j] = u1[index_mps[j][0]] / Imps // находим коэффициенты влияния для текущей точки
                }
                if (num_track_mps[i][1] == num_track_mps[j][0]) { // если номер пути конечной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1_x_find[i][j] = u2[index_mps[j][0]] / Imps // находим коэффициенты влияния для текущей точки
                }
                if (num_track_mps[i][0] == num_track_mps[j][1]) { // если номер пути начальной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2_x_find[i][j] = u1[index_mps[j][1]] / Imps // находим коэффициенты влияния для текущей точки
                }
                if (num_track_mps[i][1] == num_track_mps[j][1]) { // если номер пути конечной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2_x_find[i][j] = u2[index_mps[j][1]] / Imps // находим коэффициенты влияния для текущей точки
                }
                a_x_find[i][j] = a1_x_find[i][j] - a2_x_find[i][j]
            }
        }
    }

    /** Заполняет вектор правой части по заданному двумерному массиву с коэффициентом тока
     * для точек с одним узлом
     */
    private fun addValues_vector_b_1node(
        num_mesh: Int,
        vector_b: DoubleArray,
        arr_XI: Array<DoubleArray>,
        coeff_I: Double
    ): DoubleArray { //num_mesh -номер сетки, vector_b -вектор правой части, coeff_I - коэффициент на который умножается значение тока
        for (i in arr_XI.indices) {
            vector_b[meshes[num_mesh]!!.find_near_index_over_mesh(arr_XI[i][0])] += arr_XI[i][1] * coeff_I
        }
        return vector_b
    }

    /** Добавляет значения в вектор правой части по заданному двумерному массиву с коэффициентом тока
     * для точек с двумя узлами
     */
    private fun addValues_vector_b_2node(
        num_mesh: Int,
        vector_b: DoubleArray,
        arr_XI: Array<DoubleArray>,
        coeff_I: Double
    ): DoubleArray { //num_mesh -номер сетки, vector_b - вектор правой части А; arr_XI - двумерный массив в каждой строке  два значения: 1 - координата км, 2 - ток А; coeff_I - коэффициент тока
        var index2 = IntArray(2)
        var k: Double
        var I1: Double
        var I2: Double
        for (i in arr_XI.indices) {
            index2 = meshes[num_mesh]!!.find_2near_index_over_mesh(arr_XI[i][0]) // номер левого и правого узла
            k =
                (arr_XI[i][0] - meshes[num_mesh]!!.X[0]) / meshes[num_mesh]!!.dX // дробный индекс точки отнсительно номеров узлов сетки
            I2 = arr_XI[i][1] * coeff_I * (k - index2[0]) // ток в левый узел
            I1 = arr_XI[i][1] - I2 // ток в правый узел
            vector_b[index2[0]] += I1
            vector_b[index2[1]] += I2
        }
        return vector_b
    }

    /** Заполняет вектор правой части всех точек нулями
     */
    private fun zeros_vector_b() {
        var double_index: IntArray
        //двойной индекс
        var num_mesh: Int // номер сетки

        // проход по всем путям обнуление правой части ФОТ и ЭПС
        for (j in 0 until Ntrc) {
            // обнуление для точек в один узел сетки ФОТ
            num_mesh = tracks[j]!!.num_mesh
            for (i in tracks[j]!!.fot.indices) {
                tracks[j]!!.vector_b[meshes[num_mesh]!!.find_near_index_over_mesh(tracks[j]!!.fot[i][0])] = 0.0
            }
            //обнуление для точек в два узла сетки ЭПС
            for (i in tracks[j]!!.eps.indices) {
                double_index = meshes[num_mesh]!!.find_2near_index_over_mesh(tracks[j]!!.eps[i][0])
                tracks[j]!!.vector_b[double_index[0]] = 0.0
                tracks[j]!!.vector_b[double_index[1]] = 0.0
            }
        }

        // проход по всем МПС
        for (i in mps.indices) {
            tracks[num_track_mps[i][0]]!!.vector_b[index_mps[i][0]] = 0.0
            tracks[num_track_mps[i][1]]!!.vector_b[index_mps[i][1]] = 0.0
        }
    }

    /** Метод для расчёта мгновенной схемы ОТС
     * с указанием начального значения тока в МПС
     */
    fun calc_ots(init_I_poisk: DoubleArray): Boolean {
        return if (!verifi_I_poisk(init_I_poisk)) {
            false
        } else calc_I_poisk(init_I_poisk)
    }

    /** Метод для расчёта мгновенной схемы ОТС
     * без указания начального значения тока в поисковых точках
     */
    fun calc_ots(): Boolean {
        val init_I_poisk = DoubleArray(mps.size)
        return calc_I_poisk(init_I_poisk)
    }

    /**По сути рассчитывает токи в поисковых точках: МПС. После этого все граничные условия во всех точках с втекающим током определены
     * Расчёт каждого поискового тока сводится к вычислению невязки напряжения на данном элементе и по величине невязки корректируется ток элемента
     * По сути рассчитывается алгебраическая система уравнений методом Ньютона в цикле итераций
     */
    private fun calc_I_poisk(init_I_poisk: DoubleArray): Boolean { // возвращает истина если расчёт сошёлся, ложь в обратном случае
        //println("I_mps: ${Arrays.toString(init_I_poisk)} ")
        err_and_mes.reset_solver_error() //обнулим ошибки решателя
        err_and_mes.calc_completed = true
        verify_data(
            false,
            false,
            true
        ) // проверка исходных данных только координ ЭПС, т.к. остальное проверено при инициализации

        if (err_and_mes.data_error) { //проверка если данные корректны
            err_and_mes.solver_error = true
            err_and_mes.messeg_solver_error = "Расчёт невозможен. Ошибка исходных данных"
            return false
        }
        val N_mps = mps.size
        //количество МПС
        val U_find = DoubleArray(N_mps)
        // массивы напряжений на МПС, которые определяются всеми токами (известными и неизвестными)
        val resid_U_mps = DoubleArray(N_mps)
        var mean_resid: Double
        // средняя невязка
        val limit_mean_resid = computing_settings.convergence_U
        // задаём предельную невязку по достижении которой сходимость из класса computing_settings
        var damping_factor = computing_settings.initial_damping_factor
        //задаём коэффициент демпфирования текущее значение, на него умножается вычисленная по невязке напряжение корректирровка тока
        var mean_resid_pred: Double // значение невязки по напряжению на предыдущем шаге итераций
        var iter = 0
        val iter_max = computing_settings.max_number_iterat // счётчик итераций и максимальное число итераций
        var counter_not_exceeded: Boolean
        var convergence_not_achieved: Boolean // непревышение итераций, недостижение сходимости - булевые переменные которые определяют выход из цикла итераций
        eval_U_const() // рассчитываем напряжения на МПС от заданных токов ФОТ и ЭПС
        //println("I_2mps: ${Arrays.toString(init_I_poisk)} ")
        //println("a_x_find[i][j]: ${Arrays.deepToString(a_x_find)}")
        //нахождение токов в цикле итераций по невязке напряжения на МПС
        counter_not_exceeded = true
        convergence_not_achieved = true
        mean_resid = 100.0 //начальное значение средняя невязка до первой итерации
        while (counter_not_exceeded && convergence_not_achieved) {
            mean_resid_pred = mean_resid //предыдущая невязка обновление
            mean_resid = 0.0 //текущая невязка скидывается
            for (i in 0 until N_mps) {
                U_find[i] = U_const[i] //начинаем с постоянного напряжения (от заданных источников тока ФОТ ЭПС)

                for (j in 0 until N_mps) {
                    U_find[i] += init_I_poisk[j] * a_x_find[i][j] //добавляем напряжение от МПС
                }

                resid_U_mps[i] = U_find[i] - init_I_poisk[i] * mps[i][4] //невязка напряжения на МПС, уже изветны напряжения и Р1 и Р2
                init_I_poisk[i] += damping_factor * resid_U_mps[i] / (-0.5 * a_x_find[i][i] + mps[i][4]) //корректируем текущий поисковый ток пропорционально невязке по напряжению в этом элементе с учётом коэф. демпфирования

                mean_resid += Math.abs(resid_U_mps[i]) //обновляем невязку
            }
            mean_resid = mean_resid / N_mps //невязка именно средняя
            //если после первой итерации возрастает средняя невязка mean_resid по сравнению с ней же на предыдущей итерации mean_resid_pred, то коэффициент демпфирования в методе Ньютона уменьшаем в 0.7 раз
            if (iter > 0) {
                if (mean_resid > mean_resid_pred) {
                    damping_factor = damping_factor * 0.7
                }
            }
            iter += 1 //обновляем счётчик итераций
            counter_not_exceeded = iter < iter_max // обновляем булевые переменные выхода из цикла итераций
            convergence_not_achieved = mean_resid > limit_mean_resid
            //println("iter=$iter mean_resid=$mean_resid damping_factor=$damping_factor")
        }
        //debugLog("iter="+iter);
        computing_settings.current_state_solver = doubleArrayOf(
            (iter - 1).toDouble(),
            mean_resid,
            damping_factor
        ) // записываем текущее состояние решателя
        println("I_mps:, ${Arrays.toString(init_I_poisk)} ")
        I_mps = init_I_poisk // заносим токи в МПС в массивы родительского класса
        eval_node_from_all_I()
        zeros_vector_b()
        if (convergence_not_achieved) {
            err_and_mes.solver_error = true
            err_and_mes.messeg_solver_error =
                "Превышено максимальное число итераций равное " + computing_settings.max_number_iterat + " заданная сходимость (сред невязка по напряжен) равная " + computing_settings.convergence_U + " не достигнута."
        }
        return !convergence_not_achieved // возврат обратное к несходимости: истина если расчёт сошёлся
    }

    /*При найденных токах МПС рассчитывает напряжения и ток в рельсах решением М3ДЛ
     */
    private fun eval_node_from_all_I() {
        var g_rh: Double
        var g_lf: Double //условная проводимость слева и справа от узла сетки на схеме дискретизации рельсов
        var num_mesh: Int
        // номер сетки для начального и конечного пути МПС
        var N: Int // количество узлов сетки

        // добавляем в вектор правой часть токи МПС (токи ФОТ и ЭПС добавлены процедурой eval_U_const() )
        for (i in mps.indices) {
            tracks[num_track_mps[i][0]]!!.vector_b[index_mps[i][0]] -= I_mps[i]
            tracks[num_track_mps[i][1]]!!.vector_b[index_mps[i][1]] += I_mps[i]
        }

        //потенциал в рельсах в узлах сетки для каждого пути
        for (i in 0 until Ntrc) {
            tracks[i]!!.u = solve_3diag_band(tracks[i]!!.m3db, tracks[i]!!.vector_b)
        }

        //расчёт токов в рельсах и тока в земле в узлах сетки для каждого пути
        for (j in 0 until Ntrc) {
            num_mesh = tracks[j]!!.num_mesh //номер сетки для текущего пути
            N = meshes[num_mesh]!!.mesh_N // число узлов сетки для текущего пути
            //ток в рельсах и земле для первого узла
            g_lf = 1 / tracks[j]!!.Rv0 //проводимость слева от первого узла через волновое сопротивление в начале
            g_rh = -tracks[j]!!.m3db[0][0] //проводимость справа от первого узла через нижнюю диагональ первый элемент
            tracks[j]!!.i[0] =
                0.5 * ((0 - tracks[j]!!.u[0]) * g_lf + (tracks[j]!!.u[0] - tracks[j]!!.u[1]) * g_rh) //ток первого узла как полусумма токов слева и справа от него
            tracks[j]!!.i_grnd[0] =
                tracks[j]!!.u[0] * (tracks[j]!!.m3db[1][0] + tracks[j]!!.m3db[0][0]) //ток в земле для первого узла
            //ток в рельсах и земле для остальных узлов со второго до предпоследнего
            for (i in 1 until N - 1) {
                g_lf = -tracks[j]!!.m3db[0][i - 1] //проводимость слева от узла через нижнюю диагональ
                g_rh = -tracks[j]!!.m3db[2][i + 1] //проводимость справа от узла через верхнюю диагональ
                tracks[j]!!.i[i] =
                    0.5 * ((tracks[j]!!.u[i - 1] - tracks[j]!!.u[i]) * g_lf + (tracks[j]!!.u[i] - tracks[j]!!.u[i + 1]) * g_rh) //ток  узла как полусумма токов слева и справа от него
                tracks[j]!!.i_grnd[i] =
                    tracks[j]!!.i_grnd[i - 1] + tracks[j]!!.u[i] * (tracks[j]!!.m3db[1][i] + tracks[j]!!.m3db[0][i - 1] + tracks[j]!!.m3db[2][i + 1]) //ток в земле для узла
            }
            //ток в рельсах для последнего узла
            g_lf =
                -tracks[j]!!.m3db[2][N - 1] //проводимость слева от последнего узла через верхнюю диагональ последний элемент
            g_rh = 1 / tracks[j]!!.Rvn //проводимость справа от последнего узла через волновое сопротивление в конце
            tracks[j]!!.i[N - 1] =
                0.5 * ((tracks[j]!!.u[N - 2] - tracks[j]!!.u[N - 1]) * g_lf + (tracks[j]!!.u[N - 1] - 0) * g_rh) // ток в рельсе для псоледнего узла
            tracks[j]!!.i_grnd[N - 1] =
                tracks[j]!!.i_grnd[N - 2] + tracks[j]!!.u[N - 1] * (tracks[j]!!.m3db[1][N - 1] + tracks[j]!!.m3db[2][N - 1] - 1 / tracks[j]!!.Rvn) //ток в земле для последнего узла
        }
    }

    //------------------------------------------------------Геттеры и сеттеры головного класса + внутренние процедуры для них-------------------------------------------------------------------------
    // возвращает параметр распределенный по узлам сетки в привязке к координатам
    private fun return_X_param_rail( num_mesh: Int, X_rail: DoubleArray, arr_param_node: DoubleArray ): DoubleArray? { //num_mesh - номер сетки, X_rail[] - массив координат, км; arr_param_node - массив параметра рсапредленный по узлам сетки
        val N = X_rail.size
        var indexes: IntArray
        val out = DoubleArray(N)
        var proportion_node_left: Double
        //проверки при некорректности возврат null
        if (err_and_mes.data_error || !err_and_mes.calc_completed) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим пустой массив
            return null
        }
        if (X_rail[0] < meshes[num_mesh]!!.X[0] || X_rail[N - 1] > meshes[num_mesh]!!.X[meshes[num_mesh]!!.mesh_N - 1]) { // проверка если координаты расчётного массива за пределами сетки
            return null
        }
        //непосредственно распредление параметра в узлах сетки на массив координат в пределах сетки
        for (i in 0 until N) {
            indexes =
                meshes[num_mesh]!!.find_2near_index_over_mesh(X_rail[i]) // индексы левого и правого узла сетки относительно координаты X_rail[i]
            proportion_node_left =
                (meshes[num_mesh]!!.X[indexes[1]] - X_rail[i]) / meshes[num_mesh]!!.dX // доля для левого узла сетки обратно пропорциональна расстоянию до правого узла
            out[i] =
                arr_param_node[indexes[0]] * proportion_node_left + arr_param_node[indexes[1]] * (1 - proportion_node_left) // напряжение в левом узле на его долю + напряжение в правом узле на его долю
        }
        return out
    }

    //геттер возвращает напряжение рельс-земля в виде массива по узлам сетки
    fun get_U_rail(num_track: Int): DoubleArray { // num_track - номер пути
        return tracks[num_track]!!.u
    }

    //геттер возвращает ток в рельсах в виде массива по узлам сетки
    fun get_I_rail(num_track: Int): DoubleArray { // num_track - номер пути
        return tracks[num_track]!!.i
    }

    //геттер возвращает ток в земле в виде массива по узлам сетки для номеров путей заданных массивом
    // все пути для которых вычисляется ток в земле должны быть с единой сеткой
    fun get_I_grnd(nums_tracks: IntArray): DoubleArray { // nums_tracks - массив содержащий номера путей для которых суммировать ток в земле
        val num_mesh = tracks[nums_tracks[0]]!!.num_mesh
        // номер сетки для первого пути в массиве
        val N = meshes[num_mesh]!!.mesh_N // число узлов сетки
        val i_grnd_sum = DoubleArray(N)
        for (j in 0 until N) {
            for (i in nums_tracks.indices) {
                i_grnd_sum[j] += tracks[nums_tracks[i]]!!.i_grnd[j]
            }
        }
        return i_grnd_sum
    }

    //геттер возвращает ток в земле в виде массива по узлам сетки в сумме для всех путей
    fun get_I_grnd(): DoubleArray {
        val nums_tracks = IntArray(Ntrc)
        for (i in 0 until Ntrc) {
            nums_tracks[i] = i + 1
        }
        return get_I_grnd(nums_tracks)
    }

    //геттер возвращает напряжение рельс-земля в виде массива по заданным координатам точек
    fun get_U_rail( num_track: Int, X_rail: DoubleArray ): DoubleArray? { //  num_track - номер пути; X_rail - массив координат км;
        val num_mesh = tracks[num_track]!!.num_mesh // номер сетки для пути
        return return_X_param_rail(num_mesh, X_rail, get_U_rail(num_track))
    }

    //геттер возвращает ток в рельсах в виде массива по заданным координатам точек
    fun get_I_rail( num_track: Int, X_rail: DoubleArray ): DoubleArray? { // num_track - номер пути; X_rail - массив координат км;
        val num_mesh = tracks[num_track]!!.num_mesh // номер сетки для пути
        return return_X_param_rail(num_mesh, X_rail, get_I_rail(num_track))
    }

    //геттер возвращает ток в земле в виде массива для заданных путей по заданным координатам точек
    fun get_I_grnd( nums_tracks: IntArray, X_rail: DoubleArray ): DoubleArray? { // nums_tracks - массив содержащий номера путей для которых суммировать ток в земле; X_rail - массив координат км
        val num_mesh = tracks[nums_tracks[0]]!!.num_mesh // номер сетки для первого пути в массиве
        return return_X_param_rail(num_mesh, X_rail, get_I_grnd(nums_tracks))
    }

    //возвращает массив с токами в МПС
    fun get_I_poisk(): DoubleArray? {
        return if (err_and_mes.data_error || !err_and_mes.calc_completed) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим пустой массив
            null
        } else I_mps
    }

    /*проверка коорректности входного двумерного массива I_poisk[] содержащего токи МПС
     возвращает true если корректен, и в противном случае false
     */
    private fun verifi_I_poisk(I_poisk: DoubleArray): Boolean {
        if (I_poisk.size != mps.size) {
            err_and_mes.messeg_solver_error =
                "Ошибка при вводе массива токов МПС. Количество элементов массива должно совпадать с количеством соответствующих точек МПС"
            return false
        }
        return true
    }

    /*сеттер задает значения в МПС без выполнения непосредственного системы уравнений
      возвращает true при успешном выполнении процедуры, false  - в противном случае
      на вход принимает двумерный массив содержащий три строки, это массивы токов в следующем порядке: МПС, ЗАЗ1, ЗАЗ2
    */
    fun set_I_poisk_no_call(I_poisk: DoubleArray): Boolean {
        if (err_and_mes.data_error) { // проверка если ошибка исходных данных, то выводим  false. Без кооректного задания ФОТ, ЭПС1 и ЭПС2 - расчёт даже в этом случае невозможен
            err_and_mes.messeg_solver_error = "Задание тока в МПС не возможно. Ошибка исходных данных"
            return false
        }
        if (!verifi_I_poisk(I_poisk)) { //вызываем метод проверки корректности введенного двумерного массива I_poisk
            return false
        }
        //записываем массивы токов в поисковых точках из заданных на входе
        I_mps = I_poisk
        zeros_vector_b() // обнуление вектора правой части всех путей
        eval_U_const() // расчёт от постоянных источников
        eval_node_from_all_I() // расчёт всех узлов с учётом токов МПС
        return true
    }

    //расчитывает и выдаёт мгновенную мощность потерь в ОТС (все элементы ОТС) в кВт
    fun get_P_ots(): Double {
        if (err_and_mes.data_error || !err_and_mes.calc_completed) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим  -1
            return (-1).toDouble()
        }
        var P_ots = 0.0 // мощность ОТС накапливается в этой переменной
        for (i in 0 until Ntrc) { // проход по всем путям
            for (j in tracks[i]!!.fot.indices) { // проход по всем ФОТ данного пути, его ток с минусов
                P_ots += -tracks[i]!!.fot[j][1] * get_U_rail(i, doubleArrayOf(tracks[i]!!.fot[j][0]))!![0]
            }
            for (j in tracks[i]!!.eps.indices) { // проход по всем ЭПС данного пути
                P_ots += tracks[i]!!.eps[j][1] * get_U_rail(i, doubleArrayOf(tracks[i]!!.eps[j][0]))!![0]
            }
        }
        return P_ots
    }
}