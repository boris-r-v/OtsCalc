package ots.statistic

import ots.calc.Compute
import ots.calc.PV

/**
 * Имитирует движение поездов по участку пути
 *
 * @param cmp Класс расчета ОТС
 * @param eps Массив поездов на каждом пути: на первом пути расположены поезда из eps[0], на втором esp[1], etc. Если на пути нет поездов- то передать пустой имассив, может позже в нем поезда появятся
 * @param dXesp смещение поезда за один расчетный шаг
 */
class MoveImitator(
    internal val cmp: Compute,
    internal val eps: Array<Array<PV>>,
    private val Xeps: Double )
{
    /**
     * Проводит один этап моделирования смещая поезда на Xeps
     * проводит расчет при текущем положении поездов
     * сдвигает все поезда на Xeps
    */
    internal fun tic() {
        cmp.calcOts()
        for (tr in eps) {
            tr.forEachIndexed { i, train ->
                train.point += Xeps
                if (train.point > cmp.tracks[i].mesh.endX) {
                    train.point = cmp.tracks[i].mesh.startX
                }
            }
        }
    }
}