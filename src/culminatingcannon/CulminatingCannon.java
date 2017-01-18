
package culminatingcannon;

/**
 *
 * @author wipri9236
 */
public class CulminatingCannon {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // DECLARATIONS / INITIALZATIONS
        double v1=10, t=6,a=-9.8;
        
        
        System.out.println("dislacement is: "+dv1taOFdy(v1,t,a));
        System.out.println("initial vel is: "+dv1taOFv1(dv1taOFdy(v1,t,a),t,a));
    }
    
    public static double dxtOFv1x(double dx, double t){
        // DECLARE A CONSTANT SPEED VARIABLE
        double v1x;
        
        // CALCULATE FOR VELOCITY
        v1x = dx/t;
        
        // RETURN
        return v1x;
    }
    
    public static double dxv1xOFt(double dx, double v1x){
        // DECLARE A DELTA TIME VARIABLE
        double t;
        
        // CALCULATE FOR DELTA TIME
        t = dx/v1x;
        
        // RETURN
        return t;
    }
    
    public static double tv1xOFdx(double t, double v1x){
        // DECLARE A DELTA DISPLACEMENT
        double dx;
        
        // CALCULATE FOR DELTA DISPLACEMENT
        dx = t*v1x;
        
        // RETURN
        return dx;
    }
    
    
    // THIS IS SECTION 1 FOR PROJECTILE MOTION
    // IT USES THE FORMULA d = v1*t + 1/2 * a * t *t

    
    public static double v1ytaOFdy(double v1y,double t,double a){
        // DECLARE Δ DISPLACEMENT VARIABLE
        double dy;

        // CALCULATE FOR Δ DISPLACEMENT
        dy=v1y*t+0.5*a*Math.pow(t, 2);
        
        // RETURN
        return dy;
    }
    
    public static double dytaOFv1y(double dy,double t,double a){
        // DECLARE initial speed
        double v1y;

        // CALCULATE FOR initial speed
        v1y=(dy-(0.5*a*Math.pow(t, 2)))/t;
        
        // RETURN
        return v1y;
    }
    
    public static double dyv1ytOFa(double dy,double t,double v1y){
        // DECLARE ACCELERATION VARIABLE
        double a;

        // CALCULATE FOR initial speed
        a=(2*(dy-(v1y*t)))/Math.pow(t, 2);
        
        // RETURN
        return a;
    }
    
    public static double dyv1yaOFt(double dy,double a,double v1y){
        // DECLARE ACCELERATION VARIABLE
        double t1, t2;
        double trueTime = 133742069; // set as wild number to verify for testing
        
        // CALCULATE FOR initial speed
        t1=(-2*v1y+Math.sqrt(Math.pow(2*v1y, 2)-4*(a*(-2*dy))))/2*a;
        t2=(-2*v1y-Math.sqrt(Math.pow(2*v1y, 2)-4*(a*(-2*dy))))/2*a;
        
        // CHECK WHICH ONE OF THE TIMES IS VALID AND SET THE TRUE TIME TO IT
        if(t1<0){
           trueTime=t1; 
        }
        else if(t2<0){
            trueTime=t2;
        }
        
        return trueTime;
    }
    
    // THIS IS SECTION 1 FOR PROJECTILE MOTION
    // IT USES THE FORMULA d = v2*t - 1/2 * a * t *t
    public static double v2ytaOFdy(double v2y,double t,double a){
        // DECLARE Δ DISPLACEMENT VARIABLE
        double dy;

        // CALCULATE FOR Δ DISPLACEMENT
        dy = v2y*t - (a*t*t)/2;
        
        // RETURN
        return dy;
    }
    
    public static double dytaOFv2y(double dy,double t,double a){
        // DECLARE FINAL VELOCITY VARIBALE
        double v2y;

        // CALCULATE FOR FINAL VELOCITY
        v2 = (dy+(a*t*t)/2)/t;
        
        // RETURN
        return v2y;
    }
    
    public static double dyv2yaOFt(double dy,double v2y,double a){
        // DECLARE FINAL VELOCITY VARIBALE
        double t1, t2, trueTime = 0;

        // CALCULATE FOR FINAL VELOCITY
        t1 = (-v2y + Math.sqrt(v2y*v2y - 4*(-a/2)*(-dy)))/2*(-a);
        t2 = (-v2y + Math.sqrt(v2y*v2y - 4*(-a/2)*(-dy)))/2*(-a);
        
        // CHECK WHICH ONE OF THE TIMES IS VALID AND SET THE TRUE TIME TO IT
        if(t1<0){
           trueTime=t1; 
        }
        else if(t2<0){
            trueTime=t2;
        }
        
        return trueTime;
    }
    
    public static double dyv2ytOFa(double dy,double v2y,double t){
        // DECLARE ACCELERATION VARIABLE
        double a;

        // CALCULATE FOR initial speed
        a=(-2)*(dy-v2y*t)/(t*t);
        
        // RETURN
        return a;
    }
    
    //Method to solve for the final velocity given the initial velocity, acceleration, and time
    //Uses the formula v2=v1+at
    public static double v1yatOFv2y(double v1y, double a, double t){
        //Variable to be solved for
        double v2y;
        //Formula to solve for variable
        v2y = v1y+(a*t);
        //Return value
        return v2y;
    }
    
    //Method to solve for the initial velocity given the final velocity, acceleration, and time
    //Uses the formula v1=v2-at
    public static double v2yatOFv1y(double v2y, double a, double t){
        //Variable to be solved for
        double v1y;
        //Formula to solve for variable
        v1y = v2y-(a*t);
        //Return value
        return v1y;
    }
    
    //Method to solve for the acceleration given the initial velocity, final velocity, and time
    //Uses the formula a=(v2-v1)/t
    public static double v1yv2ytOFa(double v1y, double v2y, double t){
        //Variable to be solved for
        double a;
        //Formula to solve for variable
        a=(v2y-v1y)/t;
        //Return value
        return a;
    }
    
    //Method to solve for the time given the initial velocity, final velocity, and acceleration
    //Uses the formula t=(v2-v1)/a
    public static double v1yv2yaOFt(double v1y, double v2y, double a){
        //Variable to be solved for
        double t;
        //Formula to solve for variable
        t=(v2y-v1y)/a;
        //Return value
        return t;
    }
    
    //Method to solve for the final velocity given the initial velocity, acceleration, and displacement
    //Uses the formula v2= square root(v1*v1+2ad)
    public static double v1yadyOFv2y(double v1y, double a, double dy){
        //Variable to be solved for
        double v2y;
        //Formula to solve for variable
        v2y=Math.sqrt( (v1y*v1y) + (2*a*dy) );
        //Return value
        return v2y;
    }
    
    //Method to solve for the initial velocity given the final velocity, acceleration, and displacement
    //Uses the formula v1= square root(v2*v2-2ad)
    public static double v2yadyOFv1y(double v2y, double a, double dy){
        //Variable to be solved for
        double v1y;
        //Formula to solve for variable
        v1y=Math.sqrt( (v2y*v2y) - (2*a*dy) );
        //Return value
        return v1y;
    }
    
    //Method to solve for the acceleration given the initial velocity, final velocity, and displacement
    //Uses the formula a = (v2*v2-v1*v1)/2d
    public static double v1yv2ydyOFa(double v1y, double v2y, double dy){
        //Variable to be solved for
        double a;
        //Formula to solve for variable
        a=((v2y*v2y)-(v1y*v1y))/(2*dy);
        //Return value
        return a;
    }
    
    //Method to solve for the displacement given the initial velocity, final velocity, and acceleration
    //Uses the formula d = (v2*v2-v1*v1)/2a
    public static double v1yv2yaOFdy(double v1y, double v2y, double a){
        //Variable to be solved for
        double dy;
        //Formula to solve for variable
        dy=((v2y*v2y)-(v1y*v1y))/(2*a);
        //Return value
        return dy;
    }
    
    //Method to solve for the displacement given the initial velocity, final velocity, and time
    //Uses the formula d=((v2+v1)/2)*t 
    public static double v1yv2ytOFdy(double v1y, double v2y, double t){
        //Variable to be solved for
        double dy;
        //Formula to solve for variable
        dy=((v2y+v1y)/2)*t;
        //Return value
        return dy;
    }
    
    //Method to solve for the time given the initial velocity, final velocity, and displacement
    //Uses the formula t=2d/v2+v1
    public static double v1yv2ydyOFt(double v1y, double v2y, double dy){
        //Variable to be solved for
        double t;
        //Formula to solve for variable
        t=(2*dy)/(v2y+v1y);
        //Return value
        return t;
    }
    
    //Method to solve for the final velocity given the initial velocity, displacement, and time
    //Uses the formula v2=2d/t-v1
    public static double v1ydytOFv2y(double v1y, double dy, double t){
        //Variable to be solved for
        double v2y;
        //Formula to solve for variable
        v2y=((2*dy)/t)-v1y;
        //Return value
        return v2y;
    }
    
    //Method to solve for the initial velocity given the final velocity, displacement, and time
    //Uses the formula v1=2d/t-v2
    public static double v2ydytOFv1y(double v2y, double dy, double t){
        //Variable to be solved for
        double v1y;
        //Formula to solve for variable
        v1y=((2*dy)/t)-v2y;
        //Return value
        return v1y;
    }

 public static double v1xv1yWiththeta1(double v1x, double v1y){
        // DECLARATIONS
        double angle=0;
        double angleDegrees=0;
        
        // CALCULATE THE INITIAL ANGLE
        angle = Math.atan(v1y/v1x);
        angleDegrees=Math.toDegrees(angle);
        
        // RETURN THE DEGREES FOR THE INITAL ANGLE
        return angleDegrees; 
    }
    
    public static double v1xtheta1FORv1y(double v1x, double initialTheta){
        // DECLARATIONS
        double v1y;
        
        // CALCULATE THE
        v1y = v1x * Math.tan(initialTheta);
        
        // RETURN THE 
        return v1y; 
    }
    
    public static double v1ytheta1FORv1x(double v1y, double initialTheta){
        // DECLARATIONS
        double v1x;
        
        // CALCULATE THE
        v1x = v1x * Math.tan(initialTheta);
        
        // RETURN THE 
        return v1x; 
    }
    
    public static double v1xv2yOFv2(double v2y, double v1x){
        double v2;
        
        v2 = Math.sqrt((v2y*v2y)+(v1x*v1x));
        return v2;
        //Make sure to include theta2 in the final velocity
    }
    
    public static double v1xv1yOFv1(double v1x, double v1y){
        double v1;
        
        v1 = Math.sqrt((v1x*v1x)+(v1y*v1y));
        return v1;
        //Make sure to include theta1 in the final velocity

    
}
