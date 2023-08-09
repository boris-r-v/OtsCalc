package ots.statistic

import ots.calc.Compute
import ots.calc.PV

/**
 * Имитирует движение поездов по участку пути, т.е. перемещение ЭПС
 * Был необходим для оталадки класса обработки результатов расчета ОТС, т.е. класса Data
 *---------!!Непосредственно для расчета  ОТС НЕ_используется!!--------------
 *
 * @param cmp Класс расчета ОТС
 * @param eps Массив поездов на каждом пути: на первом пути расположены поезда из eps[0], на втором esp[1], etc. Если на пути нет поездов- то передать пустой имассив, может позже в нем поезда появятся
 * @param dEps смещение поезда за один расчетный шаг
 */
class MoveImitator(
    private val cmp: Compute,
    private val eps: Array<Array<PV>>,
    private val dEps: Double )
{
    /**
     * Проводит один этап моделирования смещая поезда на dEps
     * проводит расчет при текущем положении поездов
     * сдвигает все поезда на dEps
    */
    internal fun tic() {
        cmp.calcOts()
        for (tr in eps) {
            tr.forEachIndexed { i, train ->
                train.point += dEps
                if (train.point > cmp.tracks[i].mesh.endX) {
                    train.point = cmp.tracks[i].mesh.startX
                }
            }
        }
    }
}