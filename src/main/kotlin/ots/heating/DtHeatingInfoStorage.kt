package ots.heating
/**
 * Интерфейc, который используется для сохранения данных моделирования
 */
interface DtHeatingInfoStorage {
    /**
     * Сохранить новые данные
     * @param current - переданные на вход модели данные тока и времени
     * @param temp - значения температур элементов дроссель-трансформатора после моделирования current
     * @param sec - метки модельного времени когда сохранялись данные нагрева в отчет
     */
    fun append( current: DT_current, temp: DT_temp, sec: Double )
}
