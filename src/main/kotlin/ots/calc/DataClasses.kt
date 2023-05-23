package ots.calc

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
data class PV(  val point: Double,
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
    fun VerifyData.resetDataError() {
        dataError = false
        messegDataError = ""
    }
    fun VerifyData.resetSolverError() {
        solverError = false
        messegSolverError = ""
    }

/**
 * Класс содержит конфигурацию рассчета
 */
data class ComputingSettings( // инициализатор по умолчанию
    var convergenceU: Double = 0.01, // допустимая невязка величины напряжения усреднённая по всем точкам граничного условия, в которых происходит итерационный поиск величины тока втекающего в эти точки.
    var maxIterNumber: Int = 1000, // максимальное число итераций при расчёте величины тока, втекающего в граничные точки, по методу Ньютона в цикле.
    var initialDampingFactor: Double = 0.7, // начальное коэффициент демпфирования в методе Ньютона при расчёте величины тока, втекающего в граничные точки.
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

