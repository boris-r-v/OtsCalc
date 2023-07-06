package ots.postproc

import ots.calc.Compute
import ots.calc.PV

/**
 * Имитирует движение поездов по участку пути
 *
 * @param cmp Класс расчета ОТС
 * @param eps Массив поездов на каждом пути: на первом пути расположены поезда из eps[0], на втором esp[1], etc.
 * @param dXesp смещение поезда за один расчетный шаг
 */
class Imit(
    internal val cmp: Compute,
    internal val eps: Array<Array<PV>>,
    private val Xeps: Double,
    private val dhandler: DataProccess )
{
    /**
     * Проводит один этап моделирования смещая поезда на Xeps
    */
    internal fun tic() {
        cmp.calcOts()
        dhandler.add_one( cmp )
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