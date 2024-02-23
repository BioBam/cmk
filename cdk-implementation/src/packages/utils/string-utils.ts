export class StringUtils {
  public static printEvent(event: any) {
    console.log('--- Received Event ---');
    console.log(JSON.stringify(event, null, 2));
  }
}
