// ! "C:\Program Files\Java\jdk1.8.0_65\bin\javac" -g -source 1.6 -target 1.6 -bootclasspath "C:\Program Files\Java\jdk1.8.0_65\jre\lib\rt.jar" -extdirs "" D:\code\converter\MLDropTarget.java
import java.awt.dnd.*;
import java.awt.datatransfer.*;
import java.util.*;
import java.io.File;
import java.io.IOException;

// import java.awt.KeyEventDispatcher;
// import java.awt.KeyboardFocusManager;
// import java.awt.event.KeyEvent;

// import java.awt.Robot; 
// import java.awt.event.*;
 
public class MLDropTarget extends DropTarget
{
    /**
	 * Modified DropTarget to be used for drag & drop in MATLAB UI control.
	 */
	private static final long serialVersionUID = 1L;
//     private static boolean firstTime = true;
//     private boolean ctlDown;
    private int droptype;
	private Transferable t;
    private String[] transferData;
    
    public static final int DROPERROR = 0;
    public static final int DROPTEXTTYPE = 1;
    public static final int DROPFILETYPE = 2;
    
    @SuppressWarnings("unchecked")
    @Override
    public synchronized void drop(DropTargetDropEvent evt) {
//         ctlDown = isShiftDown(evt.getSource.getComponent(0));

// This can detect key event, but fails to detect press when focus is in file browser during drag        
//         if (firstTime) {
//             firstTime = false;
//         KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventDispatcher(new KeyEventDispatcher() {
//             
//             public boolean dispatchKeyEvent(KeyEvent ke) {
//                 switch (ke.getID()) {
//                     case KeyEvent.KEY_PRESSED:
//                         if (ke.getKeyCode() == KeyEvent.VK_CONTROL) {
//                             ctlDown = true;
//                             System.out.println(ke.getModifiers());
//                         }
//                         break;
//                     case KeyEvent.KEY_RELEASED:
//                         if (ke.getKeyCode() == KeyEvent.VK_CONTROL) {
//                             ctlDown = false;
//                         }
//                         break;
//                 }
//                 return false;
//             }
//         });
//         
//         }
        
    	// Make sure drop is accepted
    	evt.acceptDrop(DnDConstants.ACTION_COPY_OR_MOVE);
    	
    	// Set droptype to zero
    	droptype = DROPERROR;        
        
        // Get transferable and analyze
        t = evt.getTransferable();
                                
        try {
            if (t.isDataFlavorSupported(DataFlavor.javaFileListFlavor)) {
            	// Interpret as list of files
//             	List<File> fileList = (ArrayList<File>) t.getTransferData(DataFlavor.javaFileListFlavor);
            	List<File> fileList = (List<File>) t.getTransferData(DataFlavor.javaFileListFlavor); //xiangrui
            	transferData = new String[fileList.size()];
            	for (int i = 0; i < fileList.size(); i++) 
            		transferData[i] = fileList.get(i).getAbsolutePath();
            	droptype = DROPFILETYPE;
            } 
            else if (t.isDataFlavorSupported(DataFlavor.stringFlavor)) {
            	// Interpret as string            	
            	transferData = new String[1];
             	transferData[0] = (String) t.getTransferData(DataFlavor.stringFlavor);
    			droptype = DROPTEXTTYPE;
            }
            	
        } catch (UnsupportedFlavorException e) {
//         	droptype = DROPERROR;
//         	super.drop(evt);        	
//             return;
        } catch (IOException e) {
//         	droptype = DROPERROR;
//         	super.drop(evt);
//             return;
        }
        
    	// Call built-in drop method (fire MATLAB Callback)       
        super.drop(evt);
    }
    
//     public boolean isShiftDown(Component c) {
//         final List<Boolean> res = new ArrayList<Boolean>();
//         final KeyListener listener = new KeyAdapter() {
//             @Override public void keyReleased(KeyEvent e) {
//                 res.add(e.isShiftDown());
//             }
//         };
//         c.addKeyListener(listener);
//         new Robot().keyRelease(KeyEvent.VK_ALT);
//         try {Thread.sleep(50);} catch (InterruptedException e) {}
//         c.removeKeyListener(listener);
//         if (res.size() > 0)
//             return res.get(0);
//     }

    public int getDropType() {
		return droptype;
    }
//     public boolean isControlDown() {
//         return ctlDown;
//     }
	public Transferable getTransferable() {
        return t;
    }
    public String[] getTransferData() {
        return transferData;
    }
}

// class ControlDown {
//     public static boolean ctlDown;
//     public static boolean isPressed() {
//         synchronized (ControlDown.class) {
//             return ctlDown;
//         }
//     }
//     
//     public static void main(String[] args) {
//         KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventDispatcher(new KeyEventDispatcher() {
//             
//             @Override
//             public boolean dispatchKeyEvent(KeyEvent ke) {
//                 synchronized (ControlDown.class) {
// //                     ctlDown = (ke.getModifiers()==KeyEvent.CTRL_MASK);
//                     switch (ke.getID()) {
//                     case KeyEvent.KEY_PRESSED:
//                         if (ke.getKeyCode() == KeyEvent.VK_CONTROL) {
//                             ctlDown = true;
//                         }
//                         break;
// 
//                     case KeyEvent.KEY_RELEASED:
//                         if (ke.getKeyCode() == KeyEvent.VK_CONTROL) {
//                             ctlDown = false;
//                         }
//                         break;
//                     }
//                     System.out.println("ctlDown = " + ctlDown);
//                     return false;
//                 }
//             }
//         });
//     }
// }
