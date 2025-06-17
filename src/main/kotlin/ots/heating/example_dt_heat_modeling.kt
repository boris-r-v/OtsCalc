package ots.heating

class Writer(name: String = "./file.csv"): DT_temp_storage{
    val writer = java.io.PrintWriter(name)
    init {
        writer.append("current.rms,current.duration,temp_C.coil,temp_C.oil,temp_C.core,temp_C.body\n")
    }
    override fun append(current: DT_current, temp: DT_temp) {
        writer.append("${current.rms},${current.duration},${temp.coil},${temp.oil},${temp.core},${temp.body}\n")
    }
    fun close(){
        writer.close()
    }

}

/**
 * Пример расчета теплового режима работы ДТ
 */
fun main(args: Array<String>) {

    val writer = Writer()
    val dt_type = сreate_dt("DT-0.6-1000")
    val dt_hc = create_heat_transfer(20.0)
    val dt_initial_temp = create_dt_initial_temp(50.0)
    val dt_model = DT(dt_initial_temp, dt_type, dt_hc)

    val once = DT_current(1000.0, 500)
    dt_model.calc(once, writer)

    val arr = arrayOf(  DT_current(20.0, 600),
                        DT_current(10.0, 600),
                        DT_current(3000.0, 90),
                        DT_current(200.0, 600),
                        DT_current(3000.0, 90),
                        DT_current(200.0, 600),
                        DT_current(3000.0, 90),
                        DT_current(200.0, 600) )


    dt_model.calc( arr, writer)
    writer.close()
}