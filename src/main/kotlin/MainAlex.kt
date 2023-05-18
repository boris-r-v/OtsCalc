import ots.calc.OTS_DC_HeterogRail_call_kt

fun main() {
    // создадим класс ОТС с указанием конструктору количества путей и сеток
    val ots1 = OTS_DC_HeterogRail_call_kt(6, 3) //создаёмэкземпляркласса

    ots1.tracks[0]!!.fot = arrayOf(doubleArrayOf(140.5, 2300.0), doubleArrayOf(160.2, 2400.0), doubleArrayOf(176.7, 3000.0) ) // по гл1  - в классе путь 0
    ots1.tracks[4]!!.fot = arrayOf(doubleArrayOf(10.2, 1400.0)) // по отх ветке 2 путь1  - в классе путь 4
    // заполним массивы ЭПС
    ots1.tracks[0]!!.eps = arrayOf(doubleArrayOf(149.0, 800.0), doubleArrayOf(171.2, 3400.0)) // по гл1 - в классе путь 0
    ots1.tracks[1]!!.eps = arrayOf(doubleArrayOf(156.7, 1900.0) ) // погл2  - вклассепуть 1
    ots1.tracks[2]!!.eps = arrayOf(doubleArrayOf(145.3, 1600.0) ) // погл3 -  вклассепуть 2
    ots1.tracks[5]!!.eps = arrayOf(doubleArrayOf(12.1, 1400.0) ) // поотхветке 2 путь2 - вклассепуть 5
    // заполним массивы МПС: точка 1 - начальняа точка подключения, точка2 - конечная точка подключения
    // элемент массива МПС (номер пути точки1, номер пути точки2, координ точки1, координ точки2, сопротивлен)
    ots1.mps = arrayOf(doubleArrayOf(0.0, 1.0, 140.5, 140.5, 0.9e-3), doubleArrayOf(1.0, 2.0, 140.5, 140.5, 0.9e-3),
        doubleArrayOf(0.0, 1.0, 151.5, 151.5, 1.5e-3), // МПС по главным путям 1-3
        doubleArrayOf(0.0, 2.0, 155.5, 155.5, 1.6e-3),
        doubleArrayOf(0.0, 1.0, 160.2, 160.2, 1.4e-3), doubleArrayOf(1.0, 2.0, 160.2, 160.2, 1.1e-3),

        doubleArrayOf(1.0, 2.0, 167.2, 167.2, 1.4e-3),
        doubleArrayOf(0.0, 1.0, 176.7, 176.7, 0.7e-3),
        doubleArrayOf(1.0, 2.0, 177.1, 177.1, 1.8e-3),
    
        doubleArrayOf(4.0, 5.0, 10.2, 10.2, 1.0e-3), // МПС по отход2
        doubleArrayOf(0.0, 3.0, 152.5, 0.0, 1.0e-5), // соединение путь гл1 (путь 0 в классе) и однопутн отход тупик (путь 3 в классе)
        doubleArrayOf(1.0, 4.0, 170.5, 0.0, 1.0e-5), // соединение путь гл2 (путь 1 в классе) и  отход2 путь1 (путь 4 в классе)
        doubleArrayOf(2.0, 5.0, 170.5, 0.0, 1.0e-5) ) // соединение путь гл3 (путь 2 в классе) и  отход2 путь2 (путь 5 в классе)
    //заземлители и сосредоточенные продольные сопротивления учитывать не будем
    // продольные сопротивления и переходные сопротивления рельс-земля каждого пути
    for ( i in 0..2) {
        ots1.tracks[i]!!.r = arrayOf( doubleArrayOf( 189.0, 0.0254 ) ) // пути гл 1-3 продольноесопротивлен мало, переходное большое - путь бесстыковой в хорошем состоянии с высокой изоляцией от земли
        ots1.tracks[i]!!.rp = arrayOf( doubleArrayOf( 179.0, 20.0) )
    }
    ots1.tracks[3]!!.r = arrayOf( doubleArrayOf(8.0, 0.028) ) //  электрифтупикпутьзвеневой, изоляцияплохая
    ots1.tracks[3]!!.rp = arrayOf( doubleArrayOf(8.0, 1.1) )
        
    for ( i in  4..5) {
        ots1.tracks[i]!!.r = arrayOf( doubleArrayOf(15.0, 0.026) ) // двухпутнаяотходящаяветкаотх2 -путьбестыковойизоляциянорм
        ots1.tracks[i]!!.rp = arrayOf( doubleArrayOf(15.0, 10.0) )
    }
    // для однопутного тупика зададим  большое волновое сопротивление в начале и в конце линии
    ots1.tracks[3]!!.Rvn = (1e6)
    ots1.tracks[3]!!.Rv0 = (1e6)
    // для путей одходящей ветки2 зададим большое волновое сопротивление в начале линии
    ots1.tracks[4]!!.Rv0 = (1e6)
    ots1.tracks[5]!!.Rv0 = (1e6)


    // сетки
    // задание границ сетки
    ots1.meshes[0]!!.X_beg = 138.0 // сетка 0 соответствует Хгл
    ots1.meshes[0]!!.X_end = 180.0
    ots1.meshes[1]!!.X_beg = 0.0 // сетка 1  соответствует Хотх1
    ots1.meshes[1]!!.X_end = 7.0
    ots1.meshes[2]!!.X_beg = 0.0 // сетка 2 соответствует Хотх2
    ots1.meshes[2]!!.X_end = 14.0
    //присваиваем каждому пути номер сетки
    ots1.tracks[0]!!.num_mesh = 0 // сетка 0 для пути 0, 1 и 2 (гл1, гл2, гл3)
    ots1.tracks[1]!!.num_mesh = 0
    ots1.tracks[2]!!.num_mesh = 0
    ots1.tracks[3]!!.num_mesh = 1 // сетка 1 для пути 3 (однопутн тупик 7 км)
    ots1.tracks[4]!!.num_mesh = 2 // сетка 2 для пути 4 и 5 (отход двухпутн ветка)
    ots1.tracks[5]!!.num_mesh = 2


    ots1.init_OTS() // инициализируем ОТС
    if (ots1.err_and_mes.get_data_error()) { // проверим если ошибка в исходных данных
        println(ots1.err_and_mes.get_messeg_data_error()) // выведем в консоль
    }
    // рассчитаемОТС
    ots1.calc_ots()
/*
// для примера выведем ток и напряжение в рельсах вдоль всего пути номер 0 (1 гл) по координатам узлов сетки
    val mesh_track0=ots1.tracks[0].num_mesh; // номер сетки для пути 0
    val doubleX_track0 =ots1.meshes[mesh_track0].get_X(), // координаты узлов сетки для пути 0
    val U_track0 = ots1.get_U_rail(0), // напряжение в пути 0
    val I_track0 = ots1.get_I_rail(0); // токвпути 0
// По аналогии можно сделать для других путей,
// можно вторым аргументам в геттерах указывать пользовательский массив координат, но в пределах сетки.

//Также при необходимости можно вывести токи в МПС
    val doubleI_mps =ots1.get_I_poisk();
// и мгновенную мощность потерь в ОТС в Вт
    val doubleP_ots=ots1.get_P_ots();
*/
/*
    for ( i in 0..5) {
        println("track${i}.U: ${Arrays.toString(ots1.get_U_rail(i))} ")
    }
    for ( i in 0..5) {
        println("track${i}.I: ${Arrays.toString(ots1.get_I_rail(i))} ")
    }
*/
    println("P: ${ots1.get_P_ots()}")
}