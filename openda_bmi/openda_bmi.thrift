
exception ModelException {
   1: string message,
}

service BMIService 
{
    void initialize(1:string value) throws (1:ModelException error),
    void update() throws (1:ModelException error),
}