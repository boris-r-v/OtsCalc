package ots_calc

//import org.kotlinmath.*
//typealias Real = Complex

/**
 * Класс с описание одного междупутного соединителя
 *
 * @param startTrack Номер первого пути междупутного соединителя
 * @param endTrack Номер второго пути междупутного соединителя
 * @param startPoint Координата начальной точки подключения, км
 * @param endPoint Координата конечной точки подключенияб км
 * @param resValue Сопротивлнеи соединителя, Ом
 */
data class Mps( val startTrack: Int,
                val endTrack: Int,
                val startPoint: Double,
                val endPoint: Double,
                val resValue: Real,
)

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
class trackOutOfMeshException(message: String) : Exception(message)


/**
 * класс для ошибок при вычислениях, вводе данных, и сообщений о них
 */
data class verify_data( // инициализатор по умолчанию - ошибок нет, расчёт не выполнен
    var data_error: Boolean = false,
    var solver_error: Boolean = false,
    var calc_completed: Boolean = false, //data_error - ошибка в данных (к примеру, неправильная длина массива), solver_error - шибка решателя (к примеру, не достигнута сходимость, исчерпано число итераций), calc_completed  - признак что расчёт выполнен хотя бы раз
    var messeg_data_error: String = "",
    var messeg_solver_error: String = "Расчёт не выполнен", // текстовое сообщение об этих ошибках
    )    // геттеры для всех свойств
    fun verify_data.reset_data_error() {
        data_error = false
        messeg_data_error = ""
    }
    fun verify_data.reset_solver_error() {
        solver_error = false
        messeg_solver_error = ""
    }
data class Computing_settings( // инициализатор по умолчанию
    var convergence_U: Double = 0.01, // допустимая невязка величины напряжения усреднённая по всем точкам граничного условия, в которых происходит итерационный поиск величины тока втекающего в эти точки.
    var max_number_iterat: Int = 1000, // максимальное число итераций при расчёте величины тока, втекающего в граничные точки, по методу Ньютона в цикле.
    var initial_damping_factor: Double = 0.7, // начальное коэффициент демпфирования в методе Ньютона при расчёте величины тока, втекающего в граничные точки.
    var current_state_solver: DoubleArray = doubleArrayOf(0.0, 0.0, 0.0), // текущее состояние решателя содержит массив из трёх чисел полученных в конце последнего расчёта: 1 - количество итераций, 2- средняя невязка по напряжению, 3 - коэффициент демпфирования
    )

/** класс для ошибок при вычислениях, вводе данных, и сообщений о них
 */
data class Errors_and_messages ( // инициализатор по умолчанию - ошибок нет, расчёт не выполнен
    var data_error: Boolean = false,
    var solver_error: Boolean = false,
    var calc_completed: Boolean = false, //data_error - ошибка в данных (к примеру, неправильная длина массива), solver_error - шибка решателя (к примеру, не достигнута сходимость, исчерпано число итераций), calc_completed  - признак что расчёт выполнен хотя бы раз
    var messeg_data_error: String = "",
    var messeg_solver_error: String = "Расчёт не выполнен", // текстовое сообщение об этих ошибках
)
    fun Errors_and_messages.reset_data_error() {
        data_error = false
        messeg_data_error = ""
    }
    fun Errors_and_messages.reset_solver_error() {
        solver_error = false
        messeg_solver_error = ""
    }

