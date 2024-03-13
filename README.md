# Houdini_Dynamic_Microsolvers
Build Houdini dynamics solver with microsolvers




它会足够简单，可以让你省去写ui，写迭代，写几何处理等过程。但它也足够困难，因为复杂的节点堆积很容易就会产生混乱。但是完成后就可以迅速完成产品应用转化


在Houdini的SOP界面新建一个双层布料的模型，然后在DOP network节点中，创建简单的vellum solver，我们调整若干参数来来更好地模拟布料行为。
进入vellum solver节点内部，我们可以看到里面microsolver组织的非常非常复杂。不过我们从0开始用不到如此完备的功能。
现在把所有的节点都删除掉，然后我们准备开始自己使用microsolver自己新建
