package ots.heating

import kotlin.math.pow
import kotlin.math.log
import kotlin.math.abs

private val MAX_BDF_STEP = 10
class MaxBdfStepExcept(message: String) : Exception(message)
class BdfStepIsTooSmall(message: String) : Exception(message)
/**
 * Класс для хранения ошибка расчета теплового баланса ДТ на каком-то расчетном шаге
 * Учитывается отдельно ошибки для корпуса, масла, сердечника и обмотки
 */
data class CError (var oil:Double,  var coil:Double, var core:Double, var body:Double)
/**
 * Класс с данными конструкции дроссель-трансформатора
 * @param coil_HC - теплоемкость материала токовой обмотки ДТ
 * @param coil_M - масса токовой обмотки ДТ
 * @param coil_R - сопротивление обмотки тяговому току при 20Гр Цельсия
 * @param coil_S - площадь поверхности теплоотдачи обмотки
 * @param core_HC - теплопроводность материала сердечника ДТ
 * @param core_M - масса сердечника ДТ
 * @param core_S - плошать теплообмена сердечника ДТ
 * @param oil_HC - теплоемкость масла ДТ
 * @param oil_M - масса масла ДТ
 * @param body_HC - теплоемкость корпуса ДТ
 * @param body_M - масса корпуса ДТ (вместе с крышкой)
 * @param body_air_S - эквивалентная площадь поверхности (с учетом вертикальных и горизонтальных сторон ДТ) участвующая в конвективном теплообмене с воздухом
 * @param body_grayness - степень черноты поверхности ДТ
 * @param body_ground_S - площадь поверхности ДТ участвующая в кондуктивном теплообмене с землей (передача тепла от ДТ в землю)
 * @param body_oil_S - плошать теплообмена с маслом
 * @param body_radiation_S - площадь поверхности ДТ с которой возможно излучение тепла в окружающую среду (обычно площадь поверхности ДТ, без учета площади основания на котором он стоит)
 * @param body_sun_S - площадь ДТ нагреваемая солнечным светом
 * @param nominal_current - номинальный ток ДТ (пока не используется)
 * @param name - имя типа ДТ
 */
data class DT_type(
    val coil_M: Double,    val coil_HC: Double,    val coil_R:Double,    val coil_S: Double,
    val core_M: Double,    val core_HC: Double,    val core_S: Double,
    val oil_M: Double,     val oil_HC: Double,
    val body_M: Double,    val body_HC: Double,    val body_oil_S: Double,    val body_air_S: Double,    val body_radiation_S: Double,    val body_ground_S: Double,    val body_sun_S: Double,    val body_grayness: Double,
    val nominal_current: Double,
    val name: String
){}

/**
 * Класс содержит данные обратного тягового тока, протекающего через ДТ
 *
 * @param rms - Действующее значение тока в ДТ на интервале времени duration, Амперы
 * @param duration - Длительность интервала времени на котором через ДТ протекал ток со значением rms, Секунды
 */
data class DT_current( val rms: Double, val duration: Int ){}

/**
 * Класс с данными по температуре элементов дроссель-трансформатора
 * @param coil - температура обмотки
 * @param oil  - температура масла
 * @param core - температура сердечника
 * @param body - температура корпуса
 */
data class DT_temp(var coil: Double, var oil: Double, var core: Double, var body: Double ){}
/**
 * Класс с данными по коэффициентам теплопередачи и другим данным необходимым для теплового расчета
 * @param h_coil_oil - коэффициент теплопередачи от обмотки маслу, пересчитывается в зависимости от температуры масла
 * @param h_oil_core - коэффициент теплопередачи от масла к сердечнику, пересчитывается в зависимости от температуры масла
 * @param h_oil_body - коэффициент теплопередачи от масла корпусу, пересчитывается в зависимости от температуры масла
 * @param h_body_air - коэффициент теплопередачи от корпуса ДТ в воздух, пересчитывается в зависимости от температуры корпуса
 * @param h_body_ground - коэффициент теплопередачи от корпуса ДТ в землю, не пересчитывается в процессе расчета, константа
 * @param sun_radiation - поток солнечного излучения
 * @param external_temp - температура окружающей среды
 * @param cloud - облачность на небе в долях единицы, где 0-это сплошная облачность, 1-ясно, без облаков
 * @param wind - скорость ветра в м/с
 */
data class DT_heat_transfer ( var h_coil_oil: Double, var h_oil_body: Double, var h_oil_core: Double,
                              var h_body_air: Double, var h_body_ground: Double,
                              var sun_radiation: Double, var external_temp: Double,
                              var cloud: Double, var wind: Double )

/**
 * Функция для пересчета коэффициентов теплопередачи в зависимости от температуры
 */
private fun updateHtParams(dtTemp: DT_temp, dtHt: DT_heat_transfer) {
        val difT = dtTemp.body - dtHt.external_temp + 0.001
        dtHt.h_body_air = 2.1 + 1.2834 * log(difT, 2.73) + 1.51 + 4.2 * dtHt.wind
        dtHt.h_coil_oil = 0.4019 * dtTemp.oil + 110.09
        dtHt.h_oil_body = dtHt.h_coil_oil
        dtHt.h_oil_core = dtHt.h_coil_oil
}

/**
 * Класс расчета теплового баланса ДТ
 * @param temp - температура элементов ДТ
 * @param type - тип ДТ
 * @param ht - значения коэффициентов теплопередачи
 */
class DT(var temp: DT_temp, val type: DT_type, var ht: DT_heat_transfer){
    internal var temp_C = DT_temp(temp.coil-273.0, temp.oil-273.0, temp.core-273.0, temp.body-273.0)
    private var current_ = 0.0
    private val Q = mutableMapOf<String, Array<Double> >("oil" to Array(2, {0.0}), "coil" to Array(2, {0.0}), "core" to Array(2, {0.0}), "body" to Array(2, {0.0}))
    private val Csb = 5.6e-8

    /**
     * Функция дифференциального уравнения теплового баланса обмотки
     * @param - температуры: обмотки, масла, сердечника, корпуса
     * @return - рассчитанную температуру обмотки
     * f1 - сколько тепла выделилось в обмотке
     * f2 - сколько тепла отдали в масло
     */
    private fun coilOde(coilTemp: Double, oilTemp: Double, coreTemp: Double, bodyTemp:Double, idx: Int = 0 ): Double
    {
        val f1 = type.coil_R*(1+0.004*(coilTemp-293))*current_.pow(2)
        val f2 = ht.h_coil_oil*type.coil_S*(coilTemp-oilTemp)
        Q["coil"]!![idx] = f1 - f2
        //println("HEAT FOR $idx COIL ${Q["coil"]?.get(idx)}")
        return (f1 - f2) / (type.coil_M*type.coil_HC)
    }
    /**
     * Функция дифференциального уравнения теплового баланса масла
     * @param - температуры: обмотки, масла, сердечника, корпуса
     * @return - рассчитанную температуру масла
     */
    private fun oilOde (coilTemp: Double, oilTemp: Double, coreTemp: Double, bodyTemp:Double, idx: Int = 0 ): Double
    {
        val f1 = (ht.h_coil_oil*type.coil_S*(coilTemp - oilTemp) )        //heat transfer from coil to oil
        val f2 = (ht.h_oil_body*type.body_oil_S*(oilTemp - bodyTemp) )    //heat transfer from oil to body
        val f3 = (ht.h_oil_core*type.core_S*(oilTemp - coreTemp) )        //heat transfer from oil to core
        Q["oil"]!![idx] = (f1 - f2 - f3)
        //println("HEAT FOR $idx OIL ${Q["oil"]?.get(idx)}")
        return ( f1 - f2 - f3 ) / ( type.oil_M*type.oil_HC )          //oil temp
    }
    /**
     * Функция дифференциального уравнения теплового баланса сердечника
     * @param - температуры: обмотки, масла, сердечника, корпуса
     * @return - рассчитанную температуру сердечника
     */
    private fun coreOde(coilTemp: Double, oilTemp: Double, coreTemp: Double, bodyTemp:Double, idx: Int = 0 ): Double
    {
        Q["core"]!![idx] = ht.h_oil_core*type.core_S*(oilTemp - coreTemp)
        /*println("HEAT $idx FOR CORE ${Q["core"]?.get(idx)}")*/
        return ( ht.h_oil_core*type.core_S*(oilTemp - coreTemp) ) / ( type.core_M*type.core_HC )
    }
    /**
     * Функция дифференциального уравнения теплового баланса корпуса
     * @param - температуры: обмотки, масла, сердечника, корпуса
     * @return -рассчитанную температуру корпуса
     */
    private fun bodyOde (coilTemp: Double, oilTemp: Double, coreTemp: Double, bodyTemp:Double, idx: Int = 0 ): Double
    {
        val f1 = ht.h_body_air*type.body_air_S*(bodyTemp - ht.external_temp)        //конвекция с корпуса
        val f2 = ht.h_body_ground*type.body_ground_S*(bodyTemp - ht.external_temp )  //кондукция с корпуса
        val f3 = type.body_grayness*type.body_radiation_S*Csb*(bodyTemp.pow(4) - ht.external_temp.pow(4 ) ) //излучение с корпуса
        val cooling = f1 + f2 + f3
        val heating = (ht.h_oil_body*type.body_oil_S*(oilTemp - bodyTemp)) + (type.body_grayness*type.body_sun_S*ht.sun_radiation*ht.cloud)
        Q["body"]!![idx] = heating - cooling
        //println("HEAT FOR $idx BODY ${Q["body"]?.get(idx)}")
        return ( heating - cooling ) / ( type.body_M*type.body_HC )
    }
    /**
      *Расчет теплового баланса своей реализацией многоходового алгоритма численного моделирования
      * @param iError - место куда записать значения ошибки на данном шаге расчета
      * @param step - шаг расчета
      * @param maxError - разница температур между предыдущим и текущим расчетным шагом, если ошибка меньше данного числа считаем что расчет сошелся
     */
    private fun bdfIntegration(iError: CError, step: Double, maxError: Double = 1E-10 )
    {
        //println("Trace bdfIntegration**************************************")
        val error = CError(100.0,100.0,100.0,100.0)
        val dt_temp = temp
        /*Рассчитаем температуру на i+1 шаге*/
        val T_coil = dt_temp.coil + step * coilOde( dt_temp.coil, dt_temp.oil, dt_temp.core, dt_temp.body )
        val T_oil  = dt_temp.oil  + step * oilOde(  dt_temp.coil, dt_temp.oil, dt_temp.core, dt_temp.body )
        val T_core = dt_temp.core + step * coreOde( dt_temp.coil, dt_temp.oil, dt_temp.core, dt_temp.body )
        val T_body = dt_temp.body + step * bodyOde( dt_temp.coil, dt_temp.oil, dt_temp.core, dt_temp.body )

        coilOde( T_coil, T_oil, T_core, T_body,1 )
        oilOde(  T_coil, T_oil, T_core, T_body,1 )
        coreOde( T_coil, T_oil, T_core, T_body,1 )
        bodyOde( T_coil, T_oil, T_core, T_body,1 )

        var T_coil3 = dt_temp.coil + step * ((Q["coil"]!!.get(0) + Q["coil"]!!.get(1)) / 2) / (type.coil_M * type.coil_HC)
        var T_oil3 =  dt_temp.oil  + step * ((Q["oil"]!!.get(0) + Q["oil"]!!.get(1)) / 2) / (type.oil_M * type.oil_HC)
        var T_core3 = dt_temp.core + step * ((Q["core"]!!.get(0) + Q["core"]!!.get(1)) / 2) / (type.core_M * type.core_HC)
        var T_body3 = dt_temp.body + step * ((Q["body"]!!.get(0) + Q["body"]!!.get(1)) / 2) / (type.body_M * type.body_HC)

        var i = 0
        while (i < MAX_BDF_STEP && error.body > maxError && error.oil > maxError && error.coil > maxError) {

            val T_coil4 = T_coil3
            val T_oil4 = T_oil3
            val T_core4 = T_core3
            val T_body4 = T_body3

            coilOde( T_coil3, T_oil3, T_core3, T_body3,1 )
            oilOde(  T_coil3, T_oil3, T_core3, T_body3,1 )
            coreOde( T_coil3, T_oil3, T_core3, T_body3,1 )
            bodyOde( T_coil3, T_oil3, T_core3, T_body3,1 )

            T_coil3 = dt_temp.coil + step * ((Q["coil"]!!.get(0) + Q["coil"]!!.get(1)) / 2) / (type.coil_M * type.coil_HC)
            T_oil3 =  dt_temp.oil  + step * ((Q["oil"]!!.get(0) + Q["oil"]!!.get(1)) / 2) / (type.oil_M * type.oil_HC)
            T_core3 = dt_temp.core + step * ((Q["core"]!!.get(0) + Q["core"]!!.get(1)) / 2) / (type.core_M * type.core_HC)
            T_body3 = dt_temp.body + step * ((Q["body"]!!.get(0) + Q["body"]!!.get(1)) / 2) / (type.body_M * type.body_HC)
            error.coil = abs(T_coil4 - T_coil3 )
            error.oil =  abs(T_oil4 - T_oil3)
            error.core = abs(T_core4 - T_core3)
            error.body = abs(T_body4 - T_body3)
            //println("count $i, ERR: $error")
            ++i
        }
        if (i == MAX_BDF_STEP){
            throw MaxBdfStepExcept("Solution not convergent")
        }

        temp.coil = T_coil3
        temp.oil = T_oil3
        temp.core = T_core3
        temp.body = T_body3

        temp_C.coil = temp.coil-273.0
        temp_C.oil  = temp.oil-273.0
        temp_C.core = temp.core-273.0
        temp_C.body = temp.body-273.0

        iError.coil = error.coil
        iError.oil  = error.oil
        iError.core = error.core
        iError.body = error.body
    }

    /**
     * Функция расчета проводящая решение системы дифференциальных уравнений теплового баланса методом BDF
     * @param error - значение ошибки на данном шаге расчета
     * @param step - временной шаг данного расчета в секундах
     */
    private fun calcNextStep(error: CError, step: Double)
    {
        //println("calc_next_step $step")
        bdfIntegration(error, step)
        updateHtParams(temp, ht)
        //println ( "calc_next_bdf ERR : $error" )
        //println ( "calc_next_bdf Temp: $temp_C" )
        //println ( "calc_next_bdf HT  : $ht")
    }

    /**
     * Функция расчета теплового баланса на вход одно значение тока и длительности его протекания
     * @param tok - значение тока на данном временном промежутке
     * @param writer - наследник интерфейса DtHeatingInfoStorage, который умеет сохранять рассчитанные данные
     * @param writeDataEverySec - запись в отчете формируется через указанный интервал, сек (по-умолчанию 10 секунд)
     * @since Алгоритм управления шагом интегрирования:
     *  начинаем с step (1секунда) секунд, постоянно увеличивая ее в stepMult(раз)
     *  как только словили exception MaxBdfStepExcept - уменьшаем в stepMult(раз), повторяем расчет и
     *      запрещаем дальнейшее увеличение шага интегрирования - продолжаем инегрировать с этим шагом
     *  если после успешного расчета шаг окажется меньше 1 секунды и будет стоять запрет на его изменение -
     *      то увеличение шага будет разрешено
     *      если попадаем в цикл разрешаем/запрещаем - то есть ограничитель - такое можно сделать 10 раз за расчет и
     *      как только сделаем больше разрешать увеличивать шаг не будем - только можно будет уменьшать
     *  как только шаг интегрирования окажется меньше чем 1E-10, то кидаем исключение BdfStepIsTooSmall
     *      и прерываем расчет, т.к. что-то сильно пошло не так и надо разбираться с алгоритмом
     */
    fun calc( tok: DT_current, writer: DtHeatingInfoStorage, writeDataEverySec: Int = 10 )
    {
        /**
         * Описание переменных этого метода
         *  prev_write_sec - когда была сделана предыдущая запись, по времени моделирования
         *  tbeg, tend - метки времени начала и конца расчета для контроля времени выполнения
         *  step - шаг интегрирования, сек
         *  stepMult - множитель изменения шага интегрирования если получаем требуемую ошибку 1E-10, с первого захода
         *      то увеличиваем шаг на этот множитель, если нет - то уменьшаем на него же
         *  allowUpStep - флаг того что шаг можно увеличить, т.к. получили желаемую ошибку с первого раза
         *  restoredAllowUpStepTimes - сколько раз повторно разрешали увеличение шала расчета, когда он меньше 1 секунды
         *  restoredAllowUpStepMax - максимальное количество восстановлений, если превышаем то больше не разрешаем
         */
        val write_sec = 10
        var prev_write_sec = 0.0
        /*Следующий комментарий*/
        //val tbeg = System.currentTimeMillis()
        var time = 0.0
        var step = 1.0
        val stepMult = 1.2
        var allowUpStep = true
        val restoredAllowUpStepMax = 10
        var restoredAllowUpStepTimes = 0
        val error = CError(0.0,0.0,0.0,0.0)
        var iter = 0

        current_ = tok.rms
        while ( tok.duration > time )
        {
            if (step < 1E-10){
                throw BdfStepIsTooSmall(message = "DtHeating integration step is too small step<$step>")
            }
            try {
                calcNextStep(error, step)
                time += step
                //println("study time<$time>sec, step<$step>sec, temp_C $temp_C ")
            }
            catch (e: MaxBdfStepExcept)
            {
                time -= step
                step = step / stepMult
                allowUpStep = false
                //println("WARNING: ${e.message}, new step: ${step}sec ")
                continue
            }
            if ( time >= prev_write_sec + writeDataEverySec )
            {
                prev_write_sec = time
                writer.append( tok,temp_C, time)
            }
            if ( step < 1.0 && !allowUpStep && ++restoredAllowUpStepTimes < restoredAllowUpStepMax ){
                allowUpStep = true
            }
            if (allowUpStep)
            {
                step *= stepMult
            }
            iter += 1
        }

        //val tend = System.currentTimeMillis()
        //println ("total steps: $iter calc duration: ${tend-tbeg} ms" )
    }

    /**
     * Функция расчета теплового баланса на вход получает массив токов,
     * @param items - значение токов
     * @param writer - наследник интерфейса DtHeatingInfoStorage, который умеет сохранять рассчитанные данные
     * @param writeDataEverySec - запись в отчете формируется через указанный интервал, сек (по-умолчанию 10 секунд)
     *
     * @see calc
     */
    fun calc( items: Array<DT_current>, writer: DtHeatingInfoStorage, writeDataEverySec: Int = 10 ){
        items.forEach { calc(it, writer, writeDataEverySec) }
    }
    /**
     * Функция расчета теплового баланса на вход получает коллекцию токов,
     * @param items - значение токов
     * @param writer - наследник интерфейса DtHeatingInfoStorage, который умеет сохранять рассчитанные данные
     * @param writeDataEverySec - запись в отчете формируется через указанный интервал, сек (по-умолчанию 10 секунд)
     *
     * @see calc
     */
    fun calc(items: Collection<DT_current>, writer: DtHeatingInfoStorage, writeDataEverySec: Int = 10 ){
        items.forEach { calc(it, writer, writeDataEverySec) }
    }

}


