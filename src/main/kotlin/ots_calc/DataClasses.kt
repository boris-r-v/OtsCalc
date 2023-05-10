package ots_calc

import org.kotlinmath.*
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