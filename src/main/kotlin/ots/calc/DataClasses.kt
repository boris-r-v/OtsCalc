package ots.calc
import ots.complex.*
/**
 * Класс с описание одного междупутного соединителя
 *
 * @param startTrack Первый путь междупутного соединителя
 * @param endTrack Второй путь междупутного соединителя
 * @param startPoint Координата начальной точки подключения, км
 * @param endPoint Координата конечной точки подключенияб км
 * @param resValue Сопротивлнеи соединителя, Ом
 * @property startMeshIdx Номер узла сетки точки подключения соединителя к первому пути
 * @property endMeshIdx  Номер узла сетки точки подключения соединителя к первому пути
 */
data class Mps( val startTrack: Track,
                val endTrack: Track,
                val startPoint: Double,
                val endPoint: Double,
                val resValue: Real,
){
    val startMeshIdx = startTrack.mesh.findNearIndexOverMesh(startPoint)
    val endMeshIdx = endTrack.mesh.findNearIndexOverMesh(endPoint)
}


/**
 * Зачение чего то в определенной точки пути
 *
 * Применяется для
 *  1.функция погонного сопротивления рельсов вдоль пути
 *  2.функция переходного сопротивления рельсы-земля вдоль пути
 *  3.таблица ЭПС находящихся на данном пути
 *  4.таблица ЗАЗ находящихся на данном пути
 *  5.таблица сосредоточенных точечных сопротивлений в рельсах на данном пути
 *
 * @param point значение координаты точки
 * @param value значение величины в данной координате
 */
data class PV(var point: Double,
              val value: Real,
)

/**
 * Иключение генерируемое когда ордината пути вышла за границу сетки
 */
class TrackOutOfMeshException(message: String) : Exception(message)


/**
 * класс для ошибок при вычислениях, вводе данных, и сообщений о них
 */
data class VerifyData( // инициализатор по умолчанию - ошибок нет, расчёт не выполнен
    var dataError: Boolean = false,
    var solverError: Boolean = false,
    var calcCompleted: Boolean = false, //data_error - ошибка в данных (к примеру, неправильная длина массива), solver_error - шибка решателя (к примеру, не достигнута сходимость, исчерпано число итераций), calc_completed  - признак что расчёт выполнен хотя бы раз
    var messegDataError: String = "",
    var messegSolverError: String = "Расчёт не выполнен", // текстовое сообщение об этих ошибках
    )    // геттеры для всех свойств

/**
 * Класс содержит конфигурацию рассчета
 */
data class ComputingSettings(
    val convergenceU: Double = 0.01,            // допустимая невязка величины напряжения усреднённая по всем точкам граничного условия, в которых происходит итерационный поиск величины тока втекающего в эти точки.
    val maxIterNumber: Int = 100,              // максимальное число итераций при расчёте величины тока, втекающего в граничные точки, по методу Ньютона в цикле.
    val initialDampingFactor: Double = 0.7,     // начальное коэффициент демпфирования в методе Ньютона при расчёте величины тока, втекающего в граничные точки.
    val maxRelativeIter: Int = 20,              //максимамльно число итераций в расчете сходимости наведенных междуапутных напряжений, исп в callTracksUnionMesh
    val relativeConvergence: Double = 0.001,     //предел относительная сходимость среднего нааряжения на текущей и предудущей итерации в расчете наведенных междупутных напряжений, исп в callTracksUnionMesh
    val initCoefrlU: Complex = 0.049.R,            // Начальное значение коэффициента k_rlU для наведенного напряжения с track. Этот параметр меняется с текущего до 1.0
    val stepCoefrlU: Complex = 0.05.R,            // шаг возрастания коэффициента k_rlU для наведенного напряжения с track
    var currentStateSolver: DoubleArray = doubleArrayOf(0.0, 0.0, 0.0), // текущее состояние решателя содержит массив из трёх чисел полученных в конце последнего расчёта: 1 - количество итераций, 2- средняя невязка по напряжению, 3 - коэффициент демпфирования
    )

/**
 * класс для ошибок при вычислениях, вводе данных, и сообщений о них
 */
data class ErrorsAndMessages ( // инициализатор по умолчанию - ошибок нет, расчёт не выполнен
    var dataError: Boolean = false,
    var solverError: Boolean = false,
    var calcCompleted: Boolean = false, //data_error - ошибка в данных (к примеру, неправильная длина массива), solver_error - шибка решателя (к примеру, не достигнута сходимость, исчерпано число итераций), calc_completed  - признак что расчёт выполнен хотя бы раз
    var messegDataError: String = "",
    var messegSolverError: String = "Расчёт не выполнен", // текстовое сообщение об этих ошибках
)
    fun ErrorsAndMessages.resetDataError() {
        dataError = false
        messegDataError = ""
    }
    fun ErrorsAndMessages.resetSolverError() {
        solverError = false
        messegSolverError = ""
    }


/**
 * Класс ключа ассоциативного массива хранящего взаимные сопротивления между путями
 * Коммутативность достигается сохранением одного и того же массива сопротивления под
 * двумя различными колючями Key(a,b) и Key(b,a)
 * MRR - сокращение от Mesh Relative Resist
 */
data class MRRKey(val a: Track, val b: Track)
fun MRRKey.swap(): MRRKey {return MRRKey( b, a )}

/**
 * Класс хранящий набор взаимных межэдупутных споротивлений распреденных по сетке, для одной расчетной сетки
 * @param iMap ассоциативный массив междупутных сопротивьлений в исходных координатах
 * @property mesh секта к которой принадледжат все пути этого объета
 * @property data ассоциативный массив взатиных сопротивлений распределенный по узлам сетки
 */
class MeshRelativeResist(
    iMap: MutableMap<MRRKey, Array<PV>>,
)
{
    internal val mesh: Mesh = iMap.keys.first().a.mesh
    internal val data: MutableMap<MRRKey, Array<Real>> = mutableMapOf()
    init {
        iMap.forEach {
            if ( it.key.a.mesh != it.key.b.mesh ){
                throw Exception("Пути ${it.key.a.name} и ${it.key.b.name} принадлежат разным сеткам, поэтому между ними нельзя задать взаимные сопротивления, исходные данные не верны")
            }
            val dist = mesh.distributeFunctionOverMesh(it.value)
            data[it.key] = dist
            data[it.key.swap()] = dist
        }
    }

}

/**
 * Класс сожержащий массив взаимных сопротивлений путей, используется в расчетном классе
 */
class RelativeResist {
    private val arr: MutableMap<Mesh, MeshRelativeResist> = mutableMapOf()
    /**
     * Сохранить новый массив междпутуных сопротивлений
     * @param trr объект хранящий массивы междурутнеых соединений
     */
    fun set(trr: MeshRelativeResist) {
        arr[trr.mesh] = trr
    }

    /**
     * Найти и вернуть подходящий массив междупутных сопротивлений или массив нулей распреденный по узлам сетке
     * @param mesh сетка к которой принадлежат оба эти пути
     * @param tr1 первый путь
     * @param tr2 второй путь
     *
     * собственно чтобы спрятать сюда эту логику все этот класс и был создан
     * Прячем и тихо возращает массив нулей только из-за того что если он не задан то считаем что считаем без учета вляний
     * и оператор достаточно умен чтобы не забыть занести массив влияний если расчет нужен с ним
     */
    fun get(mesh: Mesh, tr1: Track, tr2: Track): Array<Real> {
        val trr = arr[mesh]
            ?: return mesh.zero     //Если для этой сетки заданы взаимные сопротивления - то вренем массив нулей
        return trr.data[MRRKey(tr1, tr2)]
            ?: throw Exception("Не зананы взаимные сопротивления между путями: ${tr1.name} и ${tr2.name}. исходные данные не верны")
    }
}