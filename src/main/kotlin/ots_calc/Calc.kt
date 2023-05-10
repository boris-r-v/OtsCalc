package ots_calc

import org.kotlinmath.R

/**
 * Класс производящий расчет токов и напряжений
 *
 * @param tracks - массив путей
 * @param meshes - массив сеток
 * @property
 */
class Calc(
    val tracks: Array<Track>,
    val meshes: Array<Mesh>,
    val mpss: Array<Mps>,
)
{
    private lateinit var I_mps: DoubleArray // массив токов междупутных соединителей МПС и токов отходящих ветвей в местах соединения с главными путями
    private lateinit var num_track_mps: Array<IntArray> // двумерный номеров путей для МПС (в каждой строке начальная и конечная точка)
    private lateinit var index_mps: Array<IntArray> //двумерный массивы  индексов узлов по сетке МПС (в каждой строке начальная и конечная точка)
    private lateinit var a_x_find: Array<Array<Real>> // матрица коэффициентов влияния тока во всех МПС на напряжения во всех МПС Ом. По главной диагонали сами на себя.
    private lateinit var U_const: Array<Real> //  массивы разности напряжений от заданных токов в МПС от начальнйо до конечной точки подключения
    init{
        init_find_index_mps()
        init_eval_a_x_find()

    }
    /**
     * Фунция проводит набор проверок корректности исходных данных
     */
    fun check(){
        println("CHECKING: ${a_x_find[0][0].toString()} ")
    }

    fun init_find_index_mps(){
        num_track_mps = Array(mpss.size) { IntArray(2) } // инициализируем массивы
        index_mps = Array(mpss.size) { IntArray(2) }
        a_x_find =Array(mpss.size) { Array<Real>(mpss.size){0.R} } // инициализируем общий массив влияния
        for (i in 0 until mpss.size) { //заполняем массивы
            num_track_mps[i][0] = mpss[i].startTrack
            num_track_mps[i][1] = mpss[i].endTrack
            index_mps[i][0] = tracks[mpss[i].startTrack].mesh.find_near_index_over_mesh(mpss[i].startPoint)
            index_mps[i][1] = tracks[mpss[i].endTrack].mesh.find_near_index_over_mesh(mpss[i].endPoint)
        }

    }
    /** Метод для решения СЛАУ с 3диагонал ленточной матрицей
     * методом двойной проходки
     * на входе матрица и вектор правой части
     * на выходе вектор ответов
     */
    fun solve_3diag_band(
        matrix_band: Array<Array<Real>>,
        vector_b: Array<Real>
    ): Array<Real> { // matrix_band – трёхдиагональная ленточная матрица, vector_b - вектор правой части
        val N = vector_b.size
        val v = Array<Real>(N){0.R}
        val u = Array<Real>(N){0.R}
        val out = Array<Real>(N){0.R}
        //прямая проходка
        v[0] = matrix_band[2][1] / -matrix_band[1][0]
        u[0] = -vector_b[0] / -matrix_band[1][0]
        for (i in 1 until N - 1) {
            v[i] = matrix_band[2][i + 1] / (-matrix_band[1][i] - matrix_band[0][i - 1] * v[i - 1])
            u[i] =
                (matrix_band[0][i - 1] * u[i - 1] - vector_b[i]) / (-matrix_band[1][i] - matrix_band[0][i - 1] * v[i - 1])
        }
        v[N - 1] = 0.0.R
        u[N - 1] =
            (matrix_band[0][N - 2] * u[N - 2] - vector_b[N - 1]) / (-matrix_band[1][N - 1] - matrix_band[0][N - 2] * v[N - 2])
        //обратная проходка
        out[N - 1] = u[N - 1]
        for (i in N - 1 downTo 1) {
            out[i - 1] = out[i] * v[i - 1] + u[i - 1]
        }
        return out
    }
    /**
     * Фукция расчета коэффициентов влиияния МПС от них самих
     */
    fun init_eval_a_x_find(){
        val N_find = this.mpss.size /* количестов поисковых точек */
        val Imps = 1000.0.R
        // уловный ток МПС для определения коэффициентов влияния
        var u1: Array<Real>
        var u2: Array<Real>
        // массивы напряжений в пути начальной и конечной точки подключения создаваемые током МПС
        val a1_x_find = Array(N_find) { Array<Real>(N_find){0.R} }
        // временные массивы влияния в МПС в начальной и
        val a2_x_find = Array(N_find) { Array<Real>(N_find){0.R} } // в конечной
        // Перебераем поисковые точки МПС
        for (i in 0 until N_find) {
            tracks[num_track_mps[i][0]].vectorB[index_mps[i][0]] = -Imps // задаём ток в начальной точке подключения в данном пути
            tracks[num_track_mps[i][1]].vectorB[index_mps[i][1]] = Imps // задаём ток в конечной точке подключения в данном пути

            u1 = solve_3diag_band(
                tracks[num_track_mps[i][0]].m3db,
                tracks[num_track_mps[i][0]].vectorB
            ) //снимаем напряжение по сетке в пути для начальной точки подключения
            u2 = solve_3diag_band(
                tracks[num_track_mps[i][1]].m3db,
                tracks[num_track_mps[i][1]].vectorB
            ) //снимаем напряжение по сетке в пути для конечной точки подключения
            tracks[num_track_mps[i][0]].vectorB[index_mps[i][0]] = 0.R // обнуляем ток в текущей начальной точке подключения
            tracks[num_track_mps[i][1]].vectorB[index_mps[i][1]] = 0.R // обнуляем ток в текущей конечной точке подключения
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
    /* Расчет постоянной составляющей напряжения в точках МПС
         */
    fun init_eval_U_const() {
        // Рассчитаем напряжение в узлах каждого пути от постоянных источников
        for (i in 0 until meshes.size) {
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
        }

        // по параметрам МПС (номерам путей и координатам точек начала и конца) запишем в U1_const и U2_const  в узлах с МПС  напряжения в рельсах от постоянных источников
        for (j in mps.indices) {
            U_const[j] =
                tracks[num_track_mps[j][0]]!!.u[index_mps[j][0]] - tracks[num_track_mps[j][1]]!!.u[index_mps[j][1]]
        }
    }

}