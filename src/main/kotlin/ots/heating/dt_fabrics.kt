package ots.heating

/**
 * Фабрика настройки тепловой модели
 * @param _ext_temp - температура внешней среды, градус Цельсия
 * @param _cloud - облачность, доли единицы, 0 - сплошная облачность, 1 - ясно, (значение по-умолчанию 0.5)
 * @param _wind - скорость ветра, м/с, (значение по-умолчанию 1.2)
 * @param _sun_radiation - удельный поток солнечной радиации на поверхность, Вт/м2 (значение по-умолчанию 800)
 */
fun create_heat_transfer(_ext_temp: Double, _cloud: Double = 0.5, _wind: Double = 1.2, _sun_radiation: Double = 800.0): DT_heat_transfer {
    val kelvin_temp = _ext_temp + 273.00001
    return DT_heat_transfer(310.0,310.0,310.0,5.35,8.69,_sun_radiation,kelvin_temp,_cloud,_wind)
}

/**
 * Фабрика создания начальной температуры ДТ
 * @param _temp - начальная температура дроссель-трансформатора, градус Цельсия
 */
fun create_dt_initial_temp( _temp: Double): DT_temp {
    val kelvin_temp = _temp + 273.005
    return DT_temp(kelvin_temp, kelvin_temp, kelvin_temp,kelvin_temp)
}
/**
 * Фабрика создания ДТ
 *  Это метод должен в боевой версии заменен на метод извлечения параметров ДТ из базы данных.
 *  Сложность в том что наполнение баы данных требует доступ к чертежам ДТ которые в открытом доступе найти не получилось
 */
fun сreate_dt(dt_type: String): DT_type
{
    if (0 == dt_type.compareTo("DT-0.6-1000") ){
        return  DT_type(40.0,390.0,0.0011,0.58,70.0, 480.0, 0.283,24.3, 1670.0, 47.0, 480.0, 0.76,0.3,1.3,0.86, 0.258,0.9,2000.0, "DT-0.6-1000")
    }
    else{
        throw Exception("Not supported type $dt_type")
    }
}